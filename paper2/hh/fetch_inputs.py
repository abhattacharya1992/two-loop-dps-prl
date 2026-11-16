#!/usr/bin/env python3
# Download the two-loop class-B gg -> HH expressions of Banerjee, Borowka, Dhani,
# Gehrmann, Ravindran, JHEP 11 (2018) 130 [arXiv:1809.05388], ancillary files
# Gauge1Simp.in and Gauge2Simp.in, and check them against the SHA-256 sums of the
# copies used for the paper. The files are not redistributed with this package.
#
# First tries the ancillary-file URLs, then the e-print source archive.
#
#   python3 fetch_inputs.py
import hashlib, io, pathlib, sys, tarfile, urllib.request

ARXIV = '1809.05388v2'
FILES = {
    'Gauge1Simp.in': '3b2f6edede722228e8526a73b1aed1feaca6ec5bb834c5632c1995a1a444b7a3',
    'Gauge2Simp.in': '6fc90c7cc0bcb2fd8025aae8ca5279a4249519873a94168ff7373e0bc5ee7a44',
}
DEST = pathlib.Path(__file__).parent/'hh_two_loop'/'anc'


def get(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'paper2-fetch/1.0'})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read()


def ok(name, data):
    return hashlib.sha256(data).hexdigest() == FILES[name]


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    missing = [n for n in FILES if not ((DEST/n).exists() and ok(n, (DEST/n).read_bytes()))]
    if not missing:
        print('ancillary files present and verified')
        return 0
    for name in list(missing):
        try:
            data = get(f'https://arxiv.org/src/{ARXIV}/anc/{name}')
            if ok(name, data):
                (DEST/name).write_bytes(data)
                missing.remove(name)
                print(f'{name}: downloaded, checksum ok')
            else:
                print(f'{name}: downloaded but checksum differs; trying the source archive')
        except Exception as e:
            print(f'{name}: direct download failed ({e}); trying the source archive')
    if missing:
        try:
            raw = get(f'https://arxiv.org/e-print/{ARXIV}')
            with tarfile.open(fileobj=io.BytesIO(raw), mode='r:*') as tar:
                for name in list(missing):
                    member = tar.extractfile(f'anc/{name}')
                    data = member.read() if member else b''
                    if data and ok(name, data):
                        (DEST/name).write_bytes(data)
                        missing.remove(name)
                        print(f'{name}: taken from the e-print archive, checksum ok')
        except Exception as e:
            print(f'e-print download failed ({e})')
    if missing:
        print('could not obtain', ', '.join(missing), '- download them by hand from arXiv:1809.05388 '
              'and place them in', DEST)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
