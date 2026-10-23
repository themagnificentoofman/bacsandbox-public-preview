"""Read-only identity adapter, run with the registered scientific interpreter."""
import json
from pathlib import Path
import sys


def inspect(checkpoint):
    source = checkpoint['source']
    root = Path(source['executing_root']) / 'code'
    sys.path.insert(0, str(root))
    from research.spatial_nutrient import runner, study_v03
    if runner.ROOT != root or runner.code_identity()['sha256'] != source['executable_sha256']:
        raise ValueError('Scientific source identity changed')
    if runner.runtime_identity() != source['runtime']:
        raise ValueError('Scientific runtime identity changed')
    grid = study_v03.verify_grid(root.parent / 'evidence/grid', verify_control_evidence=False)
    if grid['registration_sha256'] != source['registration_sha256']:
        raise ValueError('Registered grid changed')
    expected = source['candidate_manifest_bindings']
    if [r['id'] for r in grid['candidates']] != [r['id'] for r in expected]:
        raise ValueError('Candidate order changed')
    arms = 0
    for row, binding in zip(grid['candidates'], expected, strict=True):
        manifest = runner.read_manifest(Path(binding['path']))
        if manifest['manifest_sha256'] != binding['manifest_sha256']:
            raise ValueError('Candidate manifest changed')
        arms += len(manifest['seeds']) * 2
    if arms != 320:
        raise ValueError('Registered arm count changed')
    return {'status': 'registered_identity_verified', 'arms': arms,
            'scientific_sha': source['scientific_sha'], 'runtime': source['runtime'],
            'simulation_started': False}


if __name__ == '__main__':
    print(json.dumps(inspect(json.loads(Path(sys.argv[1]).read_text())), sort_keys=True))
