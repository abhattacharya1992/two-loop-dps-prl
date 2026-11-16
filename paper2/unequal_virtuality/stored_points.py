"""Read the unrounded Stage-1 helicity CSV as analysis-only records.

The CSV retains every complex helicity amplitude used in the fits. It is
not a resumable evaluation cache and omits the twenty tensor coefficients.
Those are present in fresh run_scan.py outputs and the companion raw archive.
"""
import csv
import io
import gzip
import hashlib
import json
from pathlib import Path
import mpmath as mp
from uv_core import configuration, kinematics, cpair, provenance, LOG_HELICITIES

FIELDS = ['k','dps','helicity','one_re','one_im','cf_re','cf_im','ca_re','ca_im']


def read_csv(path):
    path = Path(path)
    manifest_path = path.with_name('stage1_provenance.json')
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('schema') != 1:
        raise ValueError('unsupported stored-data manifest')
    payload = gzip.decompress(path.read_bytes()) if path.name.endswith('.gz') else path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != manifest['csv_sha256']:
        raise ValueError('stored amplitude CSV checksum mismatch')
    config = configuration(**manifest['configuration'])
    expected = {(r['k'], r['dps']): r for r in manifest['records']}
    if len(expected) != len(manifest['records']):
        raise ValueError('duplicate point metadata')
    data = {}
    with io.StringIO(payload.decode('utf-8'), newline='') as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != FIELDS:
            raise ValueError('unexpected CSV columns')
        for row in reader:
            key = (int(row['k']), int(row['dps']))
            h = row['helicity']
            if key not in expected or h not in LOG_HELICITIES:
                raise ValueError('unexpected stored point or helicity')
            point = data.setdefault(key, {})
            if h in point:
                raise ValueError('duplicate helicity row')
            point[h] = {'1l': [row['one_re'], row['one_im']],
                        '2l_CF': [row['cf_re'], row['cf_im']],
                        '2l_CA': [row['ca_re'], row['ca_im']]}
            for pair in point[h].values():
                provenance.require_pair(pair)
    if set(data) != set(expected):
        raise ValueError('missing stored points')
    records = []
    for (k, dps), H in sorted(data.items()):
        if set(H) != set(LOG_HELICITIES):
            raise ValueError('missing stored helicities')
        with mp.workdps(max(110, dps+20)):
            kin = kinematics(config, k)
            kin = {n: cpair(v) if n == 'G' else mp.nstr(v, mp.mp.dps)
                   for n, v in kin.items()}
        records.append(dict(record_type='analysis-only-helicity-projection',
                            configuration=config, k=k, dps=dps,
                            kinematics=kin, M=H,
                            original_record=expected[(k,dps)]))
    return records
