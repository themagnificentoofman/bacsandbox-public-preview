"""Run final focused and repository discovery suites in the same interpreter.

Run from the repository with project requirements and openpyxl==3.1.5 available.
This runner records real results and all skipped test identifiers/reasons.
"""
import contextlib
import importlib.util
import json
from pathlib import Path
import platform
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
OUT = ROOT / 'research/data_readiness/results'
FOCUSED = ['bacteria.test_benchmark_table', 'bacteria.test_parameter_audit',
           'bacteria.test_comparative_sweeps', 'research.data_readiness.test_zhou_tables',
           'research.data_readiness.test_prepared_observations',
           'research.data_readiness.test_s2_s13', 'research.data_readiness.test_fetch_sources',
           'research.data_readiness.test_benchmark_mapping',
           'research.data_readiness.test_feasibility']


def run_suite(name, suite):
    path = OUT / (name + '.txt')
    with path.open('w') as stream, contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
        print('Python:', sys.version)
        print('Executable:', sys.executable)
        print('Focused modules:' if name == 'final_focused_tests' else 'Discovery:', FOCUSED if name == 'final_focused_tests' else str(ROOT))
        result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    report = {'tests_run': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors),
              'skipped': [{'test': str(test), 'reason': reason} for test, reason in result.skipped],
              'successful': result.wasSuccessful(), 'log': str(path.relative_to(ROOT))}
    print(name, json.dumps({k: v for k, v in report.items() if k != 'skipped'}), 'skipped=', len(result.skipped), flush=True)
    return report


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    report = {'python': platform.python_version(), 'executable': sys.executable,
              'focused_modules': FOCUSED,
              'baseline_sha': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
              'native_extension_available': importlib.util.find_spec('bacteria_c_kernel') is not None}
    import openpyxl, et_xmlfile
    report['analysis_dependencies'] = {m.__name__: {'version': getattr(m, '__version__', 'unknown'), 'path': m.__file__} for m in (openpyxl, et_xmlfile)}
    report['focused'] = run_suite('final_focused_tests', unittest.defaultTestLoader.loadTestsFromNames(FOCUSED))
    report['full_python_discovery'] = run_suite('final_full_python_tests', unittest.defaultTestLoader.discover(start_dir=str(ROOT), top_level_dir=str(ROOT)))
    checks = {}
    for label, command in [('pip_check', [sys.executable, '-m', 'pip', 'check']), ('diff_check', ['git', 'diff', '--check'])]:
        result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
        checks[label] = {'command': command, 'exit_code': result.returncode, 'output': result.stdout + result.stderr}
    report['checks'] = checks
    report['success'] = all(report[n]['successful'] for n in ('focused', 'full_python_discovery')) and all(c['exit_code'] == 0 for c in checks.values())
    (OUT / 'final_test_run.json').write_text(json.dumps(report, indent=2) + '\n')
    return 0 if report['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
