#!/usr/bin/env python3
"""Replay the Stage-1 benchmark from retained, unrounded helicity data.

Default: checks provenance, all fit variants, two helicities, parity and
precision control without downloading VVamp. --raw-archive additionally
recontracts the twenty tensor coefficients in the companion raw archive.
This is a regression/reproducibility check, not a statistical hypothesis test.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import mpmath as mp
from uv_core import (HERE, PARTS, LOG_HELICITIES, sources, provenance,
                     predictions, complex_pair, contract, kinematics)
from stored_points import read_csv
from analyse_scan import load_records, analyse
from run_scan import validate_record


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate_provenance(manifest, retained_records=False):
    common=manifest['original_record_common_provenance']
    if not retained_records:
        require(common['sources']==provenance.file_hashes(sources()),
                'Evaluation source bytes no longer match the retained benchmark')
    else:
        # The original fingerprints are verified below against their frozen
        # metadata. Current source bytes are not represented as the old evaluator.
        require(bool(common['sources']), 'Missing historical source hashes')
    require(manifest['input_source']['archive_sha256']==
            '1c69a6beb3be85b07225af867586f242e999337a1ccf1fb940016a5341dfc091',
            'Unexpected published amplitude input')
    for r in manifest['records']:
        meta=dict(schema=1,stage='unequal-virtuality-forward',**common,
                  parameters=dict(configuration=manifest['configuration'],k=r['k'],dps=r['dps'],
                  helicities=sorted(LOG_HELICITIES),scheme='qT',mu2='s',loop_expansion='alpha_s/(2*pi)'))
        digest=hashlib.sha256(json.dumps(meta,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
        require(digest==r['fingerprint'],'Original point fingerprint inconsistent')
        require(r['diagnostics']['distinct_crossings']==4,'Missing unequal-mass crossing')
        require(r['diagnostics']['li22_bad_order']==0,'Unvalidated special-function branch')


def summarize(csv_path, retained_records=False):
    manifest=json.loads(csv_path.with_name('stage1_provenance.json').read_text())
    validate_provenance(manifest, retained_records=retained_records)
    rows=read_csv(csv_path)
    require(len(rows)==10,'Expected nine points plus one precision repeat')
    require(sorted({r['k'] for r in rows})==list(range(4,13)),'Incomplete scan')
    records,controls=load_records(csv_path)
    require(len(controls)==1 and controls[0]['k']==10 and controls[0]['dps_low']==80
            and controls[0]['dps_high']==95,'Missing 80/95-digit control')
    require(all(mp.mpf(v)<mp.mpf('1e-20') for v in controls[0]['scaled_difference'].values()),
            'Precision-control regression')
    parity=mp.mpf(0)
    for r in rows:
        for a,b in [('++++','----'),('++--','--++')]:
            for p in PARTS:
                x,y=complex_pair(r['M'][a][p]),complex_pair(r['M'][b][p])
                parity=max(parity,abs(x-y)/max(1,abs(x)))
    require(parity<mp.mpf('1e-30'),'Parity consistency regression')
    pred=predictions(records[0]['configuration'])
    frozen=json.loads((HERE/'predictions_stage1.json').read_text())
    require(frozen['configuration']==records[0]['configuration'],'Changed preregistered benchmark')
    for key,value in [('L',pred['L']),('CF_G3',pred['2l_CF'][3]),('CF_G2',pred['2l_CF'][2]),
                      ('CA_G3',pred['2l_CA'][3]),('CA_G2',pred['2l_CA'][2])]:
        require(abs(mp.mpf(frozen['benchmark_values'][key])-value)<mp.mpf('1e-45'),
                'Changed preregistered prediction: '+key)
    out=dict(configuration=records[0]['configuration'],distinct_forward_points=9,
             precision_repeat=dict(k=10,base_dps=80,repeat_dps=95),
             precision_controls=controls,parity_scaled_difference=mp.nstr(parity,15),
             logarithm='0.5*log(s*pT2/(ma2*mb2))+i*pi/2',
             assumption='r1=r2=0 two-loop boundary, tested at this configuration, not proved generally',
             primary_model=dict(leading_degree=3,power_correction_degree=3,
                                central_k=list(range(5,13)),shifted_k=list(range(4,12))),
             primary_comparisons={},model_deviations={})
    for h in ('++++','++--'):
        result=analyse(records,h)
        primary=[r for r in result['fits'] if r['power_correction_degree']==3]
        out['primary_comparisons'][h]=dict(c1=result['c1'],c1_window_shift=result['c1_window_shift'],
                                         parts={r['part']:r['comparison'] for r in primary})
        out['model_deviations'][h]={}
        for row in result['fits']:
            key=f"P={row['power_correction_degree']} {row['part']}"
            out['model_deviations'][h][key]={j:v['absolute_deviation'] for j,v in row['comparison'].items()}
        for row in primary:
            require(row['k']==list(range(5,13)),'Wrong primary fit window')
            for power,c in row['comparison'].items():
                dev,shift=mp.mpf(c['absolute_deviation']),mp.mpf(c['window_shift'])
                norm,ls=mp.mpf(c['normalization_shift']),mp.mpf(c['least_squares_shift'])
                # Broad regression thresholds; reported accuracy is the measured value.
                require(dev<(mp.mpf('1e-7') if power=='3' else mp.mpf('1e-5')),
                        f'Coefficient regression: {h} {row["part"]} G^{power}')
                require(dev<shift+norm,'Deviation exceeds observed shifted-window and normalization variations')
                require(ls<shift+norm,'Overdetermined fit leaves the observed variation range')
    return out


def verify_raw(path,manifest,csv_path):
    require(provenance.sha256_file(path)==manifest['full_tensor_archive']['sha256'],
            'Raw archive checksum mismatch')
    raw=json.loads(gzip.open(path,'rt').read())
    projected={(r['k'],r['dps']):r for r in read_csv(csv_path)}
    metadata={(r['k'],r['dps']):r for r in manifest['records']}
    require(len(raw)==len(projected),'Raw record count mismatch')
    worst=mp.mpf(0)
    for r in raw:
        validate_record(r)
        key=(r['k'],r['dps'])
        require(key in projected,'Unexpected raw point')
        require(r['M']==projected[key]['M'],'CSV projection differs from raw amplitudes')
        original=(json.dumps(r,indent=1,allow_nan=False)+'\n').encode()
        require(hashlib.sha256(original).hexdigest()==metadata[key]['raw_sha256'],
                'Raw point checksum mismatch')
        A={p:{int(j):complex_pair(v) for j,v in cols.items()} for p,cols in r['A'].items()}
        H=contract(A,kinematics(r['configuration'],r['k']),r['M'])
        for h in r['M']:
            for p in PARTS:
                expected=complex_pair(r['M'][h][p])
                err=abs(H[h][p]-expected)/max(1,abs(expected))
                worst=max(worst,err)
                require(err<mp.mpf('1e-25'),'Tensor-to-helicity reconstruction regression')
    return dict(records=len(raw),all_twenty_tensor_coefficients_recontracted=True,
                maximum_scaled_difference=mp.nstr(worst,15),
                archive_sha256=manifest['full_tensor_archive']['sha256'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-summary',action='store_true',help='write after verifying checks; otherwise compare retained summary')
    parser.add_argument('--retained-records', action='store_true',
                        help='check historical records while retaining the original source hashes')
    parser.add_argument('--raw-archive',type=Path)
    parser.add_argument('--tensor-report',type=Path)
    args=parser.parse_args()
    mp.mp.dps=110
    csv_path=HERE/'results/stage1_amplitudes.csv.gz'
    report=summarize(csv_path, retained_records=args.retained_records)
    target=HERE/'results/stage1_summary.json'
    if args.write_summary:
        provenance.atomic_write(target,report)
    else:
        require(json.loads(target.read_text())==report,'Stored summary does not reproduce exactly')
    print('PASS: nine-point scan, precision repeat, provenance, parity, unconstrained fits and retained summary')
    for h,info in report['primary_comparisons'].items():
        for part,coeffs in info['parts'].items():
            print(h,part, {j:c['absolute_deviation'] for j,c in coeffs.items()})
    if args.raw_archive:
        manifest=json.loads((HERE/'results/stage1_provenance.json').read_text())
        tensor=verify_raw(args.raw_archive,manifest,csv_path)
        if args.tensor_report:
            provenance.atomic_write(args.tensor_report,tensor)
        print('PASS: raw tensor records',tensor)


if __name__=='__main__':
    main()
