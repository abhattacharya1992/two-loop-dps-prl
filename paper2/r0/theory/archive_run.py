"""Preserve current and previous results before an explicit from-scratch run."""
import argparse
import json
from pathlib import Path
import pickle
import shutil
import tempfile

PATTERNS = {
    'gqq': ('gqq_integrands_x*.pkl', 'gqq_laurent_x*.json', 'secdec_cache_x*.pkl', 'r0_theory_quarks.json'),
    'gg': ('gg_integrands_x*.pkl', 'gg_laurent_x*.json', 'secdec_cache_gg_x*.pkl', 'r0_theory_gluons.json'),
    'n4': ('n4_integrands_x*.pkl', 'fdh_integrands_x*.pkl', 'n4_laurent_x*.json', 'fdh_laurent_x*.json',
           'secdec_cache_n4_x*.pkl', 'secdec_cache_fdh_x*.pkl', 'fdh_check_x*.json', 'kappa1_n4.json'),
}


def archive(kind, directory=Path('.')):
    directory = Path(directory)
    stored = directory / 'stored'
    stored.mkdir(exist_ok=True)
    current = sorted({p for pattern in PATTERNS[kind] for p in directory.glob(pattern)})
    backup = None

    def history_directory():
        nonlocal backup
        if backup is None:
            history = stored / 'history'
            history.mkdir(exist_ok=True)
            backup = Path(tempfile.mkdtemp(prefix=kind + '-', dir=history))
        return backup

    for path in current:
        # Interrupted tagged runs remain evidence, but must not replace a complete
        # comparison baseline. Legacy files retain their original, unverified status.
        incomplete = False
        if '_integrands_' in path.name or '_laurent_' in path.name:
            with path.open('rb' if path.suffix == '.pkl' else 'r') as stream:
                record = pickle.load(stream) if path.suffix == '.pkl' else json.load(stream)
            incomplete = isinstance(record, dict) and record.get('complete') is False
        if incomplete:
            destination = history_directory() / ('incomplete-' + path.name)
            shutil.move(str(path), destination)
            continue
        previous = stored / path.name
        if previous.exists():
            shutil.move(str(previous), history_directory() / previous.name)
        shutil.move(str(path), stored / path.name)
    if backup is not None:
        print(f'Prior versions and incomplete checkpoints preserved in {backup}')
    print(f'Preserved {len(current)} current files in {stored}; integrands will be regenerated')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kind', choices=PATTERNS)
    archive(parser.parse_args().kind)
