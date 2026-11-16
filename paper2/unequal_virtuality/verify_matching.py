#!/usr/bin/env python3
"""Check the one-loop normalisation and finite matching for two unequal-mass benchmarks.

Reanalyse retained VVamp amplitude points without imposing the coefficients.
The mean of the equal-mass quark matching function is tested at the two
available configurations; the test does not establish the formula at arbitrary masses.
"""
import argparse
import json
import mpmath as mp
from pathlib import Path
import tarfile
import tempfile

from analyse_scan import fit, load_records
from uv_core import complex_pair, kinematics, predictions

HERE = Path(__file__).resolve().parent
ARCHIVES = (
    ("lhcpheno48_20261115.tar.gz", "stage1"),
    ("lhcpheno48_benchmark2_20261115.tar.gz", "a2_b3_s17"),
)


def quark_matching(x):
    """The known equal-mass quark splitting coefficients in CF and CA."""
    lx, ly = mp.log(x), mp.log(1-x)
    z2 = mp.zeta(2)
    return 12-6*z2+2*(lx-ly)**2, 2+z2-2*(lx**2+ly**2)


def archive_records(archive, folder):
    """Load only JSON records from the named archive directory."""
    records=[]
    with tarfile.open(archive, 'r:gz') as tar:
        for member in tar.getmembers():
            if member.isfile() and member.name.startswith(folder+'/') and member.name.endswith('.json'):
                if '/' in member.name[len(folder)+1:]:
                    raise ValueError('Unexpected nested archive path')
                records.append(json.loads(tar.extractfile(member).read()))
    if len(records)!=10:
        raise ValueError(f'Expected nine forward points and one precision repeat: found {len(records)}')
    with tempfile.TemporaryDirectory() as directory:
        for i, record in enumerate(records):
            (Path(directory)/f'uv_{i:02d}.json').write_text(json.dumps(record))
        chosen, controls = load_records(Path(directory))
    if sorted(r['k'] for r in chosen)!=list(range(4,13)) or len(controls)!=1:
        raise ValueError('Missing forward point or precision control')
    return chosen, controls


def compute():
    results=[]
    with mp.workdps(110):
        for filename, folder in ARCHIVES:
            archive=HERE/'results/local_reproductions'/filename
            records,controls=archive_records(archive, folder)
            kin=kinematics(records[0]['configuration'],5)
            a,b,s=(kin[key] for key in ('ma2','mb2','s'))
            x1,x2=kin['x1'],kin['x2']
            Lx=mp.log(mp.sqrt(a*b)/s)
            Lg=-2*Lx+1j*mp.pi
            lm1,lm2=mp.log(a/s)-1j*mp.pi,mp.log(b/s)-1j*mp.pi
            # One-loop hard-current corrections for the two different boson masses.
            h_cf=-(lm1**2+lm2**2)+3*(lm1+lm2)-16+2*mp.zeta(2)
            # Coefficients of CF, CA in the kernel contribution in the common G variable.
            kernel_cf=2*Lg**2+4*Lx*Lg+6*Lg
            kernel_ca=-Lg**2-4*Lx*Lg
            m1=quark_matching(x1);m2=quark_matching(x2)
            expected=[(m1[i]+m2[i])/2 for i in range(2)]
            predicted_logs=predictions(records[0]['configuration'])
            leading_deviation=matching_deviation=window_change=lower_model_deviation=mp.mpf(0)
            row={'configuration':records[0]['configuration'],
                 'x1':mp.nstr(x1,20),'x2':mp.nstr(x2,20),
                 'expected_CF':mp.nstr(expected[0],18),'expected_CA':mp.nstr(expected[1],18),
                 'precision_control':controls[0], 'helicities':{}}
            for h in ('++++','++--'):
                c1fit=fit(records,[complex_pair(r['M'][h]['1l']) for r in records],1,2)
                if c1fit is None:raise ValueError('No one-loop coefficient fit')
                c1=c1fit['coeff'][1]
                c1prediction=-8j*mp.pi*(x1*(1-x2) if h=='++++' else (1-x1)*x2)
                c1error=abs(c1-c1prediction)
                if c1error>mp.mpf('1e-11')*abs(c1prediction):
                    raise AssertionError('One-loop splitting coefficient does not match')
                values={p:[complex_pair(r['M'][h][p]) for r in records]
                        for p in ('2l_CF','2l_CA')}
                out={'c1_prediction':mp.nstr(c1prediction,19),
                     'c1_absolute_deviation':mp.nstr(c1error,16)}
                for part, idx in (('2l_CF',0),('2l_CA',1)):
                    f=fit(records,values[part],3,3)
                    shifted=fit(records,values[part],3,3,window=1)
                    lower=fit(records,values[part],3,2)
                    if not all(z is not None for z in (f,shifted,lower)):
                        raise ValueError('Incomplete two-loop fitting windows')
                    for j in (2,3):
                        ratio=f['coeff'][j]/c1
                        leading_deviation=max(leading_deviation, abs(ratio-predicted_logs[part][j]))
                        window_change=max(window_change, abs(ratio-shifted['coeff'][j]/c1))
                        lower_model_deviation=max(lower_model_deviation,
                                                  abs(ratio-lower['coeff'][j]/c1))
                    # Source coefficients use alpha_s/(2*pi); manuscript uses alpha_s/(4*pi).
                    c11=2*f['coeff'][1]/c1
                    c11shift=2*shifted['coeff'][1]/c1
                    c11lower=2*lower['coeff'][1]/c1
                    # The qT remainder differs from the minimally subtracted remainder
                    # by -CA*zeta_2 times the one-loop amplitude.
                    correction=(-h_cf-kernel_cf) if idx==0 else (mp.zeta(2)-kernel_ca)
                    r0=c11+correction
                    difference=abs(r0-expected[idx]); shift=abs(c11-c11shift)
                    matching_deviation=max(matching_deviation,difference)
                    model_change=abs(c11-c11lower)
                    if difference>mp.mpf('1e-8'):
                        raise AssertionError(f'{h} {part}: finite matching discrepancy {difference}')
                    out[part]={'extracted_real':mp.nstr(mp.re(r0),18),
                               'extracted_imag':mp.nstr(mp.im(r0),16),
                               'complex_deviation':mp.nstr(difference,16),
                               'fit_window_shift':mp.nstr(shift,16),
                               'change_for_P2':mp.nstr(model_change,16)}
                row['helicities'][h]=out
            row['leading_log_max_deviation']=mp.nstr(leading_deviation,17)
            row['leading_log_window_shift']=mp.nstr(window_change,17)
            row['leading_log_P2_change']=mp.nstr(lower_model_deviation,17)
            row['finite_matching_max_deviation']=mp.nstr(matching_deviation,17)
            results.append(row)
    return {'description':'Independent extraction from retained unequal-mass VVamp points; no new amplitudes',
            'interpretation':'Two numerical tests of the mean matching expression, not a general proof',
            'results':results}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=compute()
    path=HERE/'results/matching_checks.json'
    if args.write:path.write_text(json.dumps(result,indent=2)+'\n')
    if args.check:
        if json.loads(path.read_text())!=result:
            raise AssertionError('Stored unequal-mass matching results are not reproduced')
    for row in result['results']:
        for h, values in row['helicities'].items():
            print(row['configuration'],h,'c1 error',values['c1_absolute_deviation'],
                  'r0 errors',values['2l_CF']['complex_deviation'],values['2l_CA']['complex_deviation'])
    print('PASS: unequal-mass matching and one-loop normalisation')


if __name__=='__main__':main()
