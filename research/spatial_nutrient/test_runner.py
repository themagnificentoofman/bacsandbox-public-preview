"""Failure persistence and registered-denominator integration checks."""
import copy
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from . import runner
from .protocol import load_protocol


def small_protocol():
    p=load_protocol()
    p['species'][0]['init']=1
    p['species'][1]['init']=2
    p['clock'].update(horizon_model_seconds=1.0,steps=10)
    return p


def test_registration_saves_inputs_and_rejects_overwrite(tmp_path):
    out=tmp_path/'case'
    manifest=runner.register(out,small_protocol(),'diagnostic',[45])
    assert len(manifest['pairs'])==1
    assert runner.input_path(out,45,'U').name=='U.json.gz'
    assert runner.input_path(out,45,'H').name=='H.json.gz'
    assert runner.read_manifest(out)['manifest_sha256']==manifest['manifest_sha256']
    with pytest.raises(FileExistsError):runner.register(out,small_protocol(),'diagnostic',[45])


def test_failure_retained_and_retry_uses_same_saved_input(tmp_path):
    out=tmp_path/'case'
    manifest=runner.register(out,small_protocol(),'diagnostic',[46])
    with patch.object(runner,'build_simulation',side_effect=RuntimeError('injected setup fault')):
        failed=runner.execute_arm(out,46,'U')
    assert failed['status']=='failed'
    assert failed['reached_time_s']==0
    with pytest.raises(RuntimeError,match='retry required'):
        runner.execute_arm(out,46,'U')
    completed=runner.execute_arm(out,46,'U',retry=True)
    assert completed['status']=='completed',completed
    assert completed['input_sha256']==failed['input_sha256']
    assert (out/'runs/46/U/attempt-001/result.json').exists()
    assert completed['attempt']==2
    with pytest.raises(ValueError,match='missing completed arm'):
        runner.collect(out)
    assert runner.execute_arm(out,46,'H')['status']=='completed'
    _,records=runner.collect(out)
    assert len(records)==2


def test_tampering_and_main_without_freeze_rejected(tmp_path):
    out=tmp_path/'case'
    runner.register(out,small_protocol(),'diagnostic',[47])
    path=runner.input_path(out,47,'U')
    data=runner.read_inputs(out,47,'U');data['orientation_rad']=123
    with runner._compressed_jsonl(path) as handle:
        handle.write(json.dumps(data))
    with pytest.raises(ValueError,match='inputs changed'):runner.execute_arm(out,47,'U')
    with pytest.raises(ValueError,match='qualification freeze'):
        runner.register(tmp_path/'main',load_protocol(),'main')
    with pytest.raises(ValueError,match='fresh main seeds'):
        runner.register(tmp_path/'bad',small_protocol(),'diagnostic',[1000])


def test_missing_registration_cannot_run(tmp_path):
    with pytest.raises(FileNotFoundError):runner.run(tmp_path,workers=1)


def test_completed_results_validate_physical_endpoint_and_analysis(tmp_path):
    out=tmp_path/'case'
    runner.register(out,small_protocol(),'diagnostic',[48,49])
    for seed in (48,49):
        for arm in ('U','H'):
            result=runner.execute_arm(out,seed,arm)
            assert result['status']=='completed',result
            assert result['reached_time_s']==1.0
    result=runner.analyze(out)
    assert result['n_pairs']==2


def test_duplicate_uniform_inputs_replay_exactly(tmp_path):
    records=[];states=[]
    for name in ('first','second'):
        out=tmp_path/name
        runner.register(out,small_protocol(),'diagnostic',[50])
        result=runner.execute_arm(out,50,'U')
        assert result['status']=='completed'
        records.append(result)
        states.append(json.loads((out/'runs/50/U/attempt-001/rng_final.json').read_text()))
    assert records[0]==records[1]
    assert states[0]==states[1]


def test_zero_contrast_ecological_replay_is_exact(tmp_path):
    p=small_protocol();p['feeding']['contrast_nodes']=0
    out=tmp_path/'null'
    runner.register(out,p,'diagnostic',[51])
    u=runner.execute_arm(out,51,'U');h=runner.execute_arm(out,51,'H')
    assert u['status']==h['status']=='completed'
    for key in ('counts','focal_count','extinction_time_s','ledger','diagnostics','tail_min_focal_count'):
        assert u[key]==h[key]


def test_interrupted_completion_marker_repaired_and_outputs_verified(tmp_path):
    out=tmp_path/'marker'
    runner.register(out,small_protocol(),'diagnostic',[52])
    original=runner.execute_arm(out,52,'U')
    assert original['status']=='completed'
    marker=out/'runs/52/U/completed.json'
    marker.unlink()
    with pytest.raises(RuntimeError,match='completion marker'):
        runner.execute_arm(out,52,'U')
    assert runner.execute_arm(out,52,'U',retry=True)==original
    assert marker.exists()
    events=runner.artifact_path(out/'runs/52/U/attempt-001', 'events.jsonl')
    events.write_text('truncated\n')
    with pytest.raises(ValueError,match='integrity mismatch'):
        runner.execute_arm(out,52,'U',retry=True)


def test_dead_worker_lock_recovery_retains_evidence(tmp_path):
    out=tmp_path/'interrupted'
    runner.register(out,small_protocol(),'diagnostic',[53])
    arm=out/'runs/53/U';arm.mkdir(parents=True)
    runner.atomic_json(arm/'active.lock',{'pid':999999999,'host':runner.platform.node()})
    assert runner.recover_interrupted(out,53,'U')
    assert (arm/'recovery-001.json').exists()
    assert runner.execute_arm(out,53,'U',retry=True)['status']=='completed'


def test_compressed_evidence_roundtrip_and_ambiguous_copy_rejected(tmp_path):
    out=tmp_path/'compressed'
    runner.register(out,small_protocol(),'diagnostic',[54])
    record=runner.execute_arm(out,54,'U')
    assert record['status']=='completed',record
    attempt=out/'runs/54/U/attempt-001'
    events=runner.artifact_path(attempt,'events.jsonl')
    assert events.name=='events.jsonl.gz'
    with runner.open_jsonl(events) as handle:
        parsed=[json.loads(line) for line in handle]
    assert any(e['kind']=='feed' for e in parsed)
    assert sum(e.get('amount_pg',0) for e in parsed if e['kind']=='feed')==record['ledger']['feed_pg']
    (attempt/'events.jsonl').write_text('unregistered alternate copy\n')
    with pytest.raises(ValueError,match='Expected one retained artifact'):
        runner.execute_arm(out,54,'U',retry=True)


def test_historical_plain_inputs_and_output_logs_remain_readable(tmp_path):
    out = tmp_path / 'historical'
    runner.register(out, small_protocol(), 'diagnostic', [55])
    expected = []
    for arm in ('U', 'H'):
        # Recreate the historical on-disk format before execution: input identity
        # binds parsed data, while retained output identity binds the file bytes.
        compressed_input = runner.input_path(out, 55, arm)
        inputs = runner.read_inputs(out, 55, arm)
        runner.atomic_json(compressed_input.with_suffix(''), inputs)
        compressed_input.unlink()
        assert runner.read_inputs(out, 55, arm) == inputs
        record = runner.execute_arm(out, 55, arm)
        assert record['status'] == 'completed', record
        attempt = out / 'runs' / '55' / arm / 'attempt-001'
        for name in ('events.jsonl', 'trajectory.jsonl'):
            compressed_log = runner.artifact_path(attempt, name)
            with runner.open_jsonl(compressed_log) as handle:
                text = handle.read()
            plain_log = compressed_log.with_suffix('')
            plain_log.write_text(text, encoding='utf-8')
            compressed_log.unlink()
            assert runner.artifact_path(attempt, name) == plain_log
            with runner.open_jsonl(plain_log) as handle:
                assert handle.read() == text
        # Emulate a historical producer's own completion hashes; this is only
        # fixture construction, never a migration of retained research evidence.
        record['outputs'] = runner.artifact_manifest(attempt)
        runner.atomic_json(attempt / 'result.json', record)
        runner.atomic_json(attempt.parent / 'completed.json',
                           {'attempt': 1, 'result_sha256': runner.digest(record)})
        assert set(record['outputs']) == {'events.jsonl', 'trajectory.jsonl', 'rng_final.json'}
        assert runner.execute_arm(out, 55, arm) == record
        expected.append(record)
    _, records = runner.collect(out)
    assert records == expected
