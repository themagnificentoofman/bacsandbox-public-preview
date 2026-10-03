"""Fetch or verify pinned public supplements; never overwrite differing local bytes.

Bibliographic attribution and redistribution boundaries: DOCUMENTATION.md §10.6.
Network access occurs only when this CLI is explicitly run without --verify-only.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
import urllib.request

HERE = Path(__file__).resolve().parent
MANIFESTS = ['supplements', 'supplements_s2_s13']


def verify_bytes(data, entry):
    if len(data) != entry['size'] or hashlib.sha256(data).hexdigest() != entry['sha256'] or hashlib.md5(data).hexdigest() != entry['md5']:
        raise ValueError(f"source integrity failure: {Path(entry['key']).name}")


def materialize(entry, destination, verify_only=False, opener=urllib.request.urlopen):
    destination = Path(destination)
    if destination.exists():
        verify_bytes(destination.read_bytes(), entry)
        return 'verified_existing'
    if verify_only:
        raise FileNotFoundError(f'missing source: {destination}; run fetch_sources or supply the pinned original')
    with opener(entry['url'], timeout=30) as response:
        data = response.read(entry['size'] + 1)
    verify_bytes(data, entry)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=destination.parent, prefix='.source-', delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(data)
        # Exclusive installation protects an existing file, including races.
        try:
            os.link(temporary, destination)
        except FileExistsError:
            verify_bytes(destination.read_bytes(), entry)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return 'downloaded_verified'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-only', action='store_true')
    parser.add_argument('--output-dir', type=Path, default=HERE / 'raw')
    args = parser.parse_args()
    results = []
    try:
        for group in MANIFESTS:
            manifest = json.loads((HERE / 'provenance' / (group + '.json')).read_text())
            for entry in manifest['files']:
                name = Path(entry['key']).name
                state = materialize(entry, args.output_dir / group / name, args.verify_only)
                results.append({'source': name, 'status': state})
    except (OSError, ValueError) as error:
        parser.exit(2, f'{error}\n')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
