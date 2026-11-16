"""Fetch and authenticate the independent six-point amplitude ancillary data."""
from pathlib import Path
import hashlib
import io
import urllib.request
import zipfile

HERE=Path(__file__).resolve().parent
URL='https://s3df.slac.stanford.edu/people/lance/Cosmic/hexMRKL1-7.zip'

def main():
    out=HERE/'inputs'
    hashes={name:digest for digest,name in
            (line.split() for line in (out/'SHA256SUMS').read_text().splitlines())}
    with urllib.request.urlopen(URL,timeout=60) as response:
        data=response.read()
    if hashlib.sha256(data).hexdigest()!=hashes['hexMRKL1-7.zip']:
        raise ValueError('Primary archive checksum mismatch')
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        entries=[n for n in archive.namelist() if Path(n).name=='hexMRKL1-7.m']
        if len(entries)!=1:raise ValueError('Ambiguous primary amplitude file')
        primary=archive.read(entries[0])
    if hashlib.sha256(primary).hexdigest()!=hashes['hexMRKL1-7.m']:
        raise ValueError('Primary amplitude checksum mismatch')
    (out/'hexMRKL1-7.m').write_bytes(primary)
    print('Authenticated hexMRKL1-7.m')

if __name__=='__main__':main()
