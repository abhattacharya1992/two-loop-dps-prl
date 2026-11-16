#!/usr/bin/env python3
"""C++ build-cache decisions with a synthetic compiler; no VVamp calculation."""
import json
import os
import pathlib
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import check_cxx


COMPILER = r'''#!/usr/bin/env python3
import hashlib, json, pathlib, sys
here = pathlib.Path(__file__).resolve().parent
args = sys.argv[1:]
if args == ['--version']:
    print('synthetic compiler 1')
elif args == ['-dumpmachine']:
    print('synthetic-target')
elif len(args) == 1 and args[0].startswith('-print-prog-name='):
    print(__file__)
elif args == ['-print-file-name=libquadmath.so']:
    print(here / 'libquadmath.so')
else:
    action = 'compile' if '-c' in args else 'link'
    with (here / 'calls.jsonl').open('a') as stream:
        stream.write(json.dumps({'action': action, 'args': args}) + '\n')
    if (here / ('fail-' + action)).exists():
        print('synthetic ' + action + ' failure', file=sys.stderr)
        sys.exit(7)
    output = pathlib.Path(args[args.index('-o') + 1])
    if action == 'compile':
        source = pathlib.Path(args[args.index('-c') + 1])
        header = here / 'external-header.h'
        digest = hashlib.sha256(source.read_bytes() + header.read_bytes()).hexdigest()
        output.write_text('object ' + digest)
        dependency = pathlib.Path(args[args.index('-MF') + 1])
        dependency.write_text(str(output) + ': ' + str(source) + ' ' + str(header) + '\n')
    else:
        objects = [pathlib.Path(arg) for arg in args if arg.endswith('.o')]
        output.write_text('linked ' + ''.join(obj.read_text() for obj in objects))
'''


class BuildCacheTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = pathlib.Path(self.temporary.name)
        self.tools = self.root / "tools"
        self.tools.mkdir()
        self.compiler = self.tools / "compiler"
        self.compiler.write_text(COMPILER)
        self.compiler.chmod(0o755)
        (self.tools / "libquadmath.so").write_bytes(b"synthetic library")
        (self.tools / "external-header.h").write_text("header version 1")
        self.source = self.root / "work" / "src"
        self.source.mkdir(parents=True)
        self.files = ["one.cpp", "two.cpp"]
        for name in self.files:
            (self.source / name).write_text("synthetic source " + name)
        for name in ("ggvvamp.xh", "ggvvamp.xi", "calls.inc"):
            (self.source / name).write_text("synthetic " + name)
        self.compiler_patch = patch.object(check_cxx, "CXX", str(self.compiler))
        self.compiler_patch.start()

    def tearDown(self):
        self.compiler_patch.stop()
        self.temporary.cleanup()

    def build(self, **kwargs):
        return check_cxx.compile_all(self.source, self.files, **kwargs)

    def calls(self):
        path = self.tools / "calls.jsonl"
        return [json.loads(line) for line in path.read_text().splitlines()] if path.exists() else []

    def test_matching_build_is_reused_without_compilation_or_linking(self):
        exe = self.build()
        calls = self.calls()
        self.assertEqual(sum(c["action"] == "compile" for c in calls), 2)
        self.assertEqual(sum(c["action"] == "link" for c in calls), 1)
        self.assertEqual(self.build(), exe)
        self.assertEqual(self.calls(), calls)

    def test_source_header_compiler_library_and_flags_changes_trigger_rebuild(self):
        self.build()
        for path in (self.source / "one.cpp", self.source / "ggvvamp.xh",
                     self.tools / "external-header.h", self.compiler, self.tools / "libquadmath.so"):
            with self.subTest(path=path.name):
                previous_calls = len(self.calls())
                old_mtime = path.stat().st_mtime_ns
                with path.open("a") as stream:
                    stream.write("\n# changed\n")
                os.utime(path, ns=(old_mtime, old_mtime))
                self.build()
                self.assertEqual(len(self.calls()) - previous_calls, 3)
        previous_calls = len(self.calls())
        self.build(flags=["-O0", "-std=gnu++11"])
        self.assertEqual(len(self.calls()) - previous_calls, 3)

    def test_legacy_tampered_or_missing_products_are_rebuilt(self):
        for name in ("one.o", "two.o", "ggvv_cxx"):
            (self.source / name).write_text("legacy output")
        self.build()
        self.assertEqual(len(self.calls()), 3)
        (self.source / "one.o").write_text("changed object")
        self.build()
        self.assertEqual(len(self.calls()), 6)
        (self.source / "ggvv_cxx").unlink()
        self.build()
        self.assertEqual(len(self.calls()), 9)

    def test_retargeted_header_symlink_invalidates_the_build(self):
        header = self.tools / "external-header.h"
        first, second = self.tools / "header-one.h", self.tools / "header-two.h"
        first.write_text("first header")
        second.write_text("second header")
        header.unlink()
        header.symlink_to(first)
        self.build()
        before = len(self.calls())
        header.unlink()
        header.symlink_to(second)
        self.build()
        self.assertEqual(len(self.calls()) - before, 3)

    def test_compile_and_link_failures_preserve_previous_build(self):
        self.build()
        paths = [self.source / name for name in ("one.o", "two.o", "ggvv_cxx", "build_manifest.json")]
        before = {p.name: p.read_bytes() for p in paths}
        (self.source / "one.cpp").write_text("changed source")
        for phase in ("compile", "link"):
            with self.subTest(phase=phase):
                failure = self.tools / ("fail-" + phase)
                failure.touch()
                with self.assertRaises(subprocess.CalledProcessError) as caught:
                    self.build()
                self.assertEqual(caught.exception.returncode, 7)
                self.assertEqual({p.name: p.read_bytes() for p in paths}, before)
                failure.unlink()

    def test_prepare_refreshes_existing_original_sources_and_headers_by_content(self):
        original = self.root / "original"
        original.mkdir()
        # Cheap text fixtures cover the real list of generated files without compiling it.
        for k in range(3, 2464):
            (original / f"p{k:04d}.cpp").write_text("// original\n")
        (original / "p2464.cpp").write_text("T CF = n(4,3); T CA = n(3);\n")
        for name in ("ggvvamp.xh", "ggvvamp.xi"):
            (original / name).write_text("// original header\n")
        copied, _ = check_cxx.prepare(original, self.root / "prepared")
        for name in ("p0003.cpp", "ggvvamp.xh", "ggvvamp.xi"):
            (original / name).write_text("// updated\n")
        check_cxx.prepare(original, self.root / "prepared")
        for name in ("p0003.cpp", "ggvvamp.xh", "ggvvamp.xi"):
            self.assertEqual((copied / name).read_bytes(), (original / name).read_bytes())
        self.assertIn("harness_CF", (copied / "p2464.cpp").read_text())

    def test_make_cxx_failure_returns_nonzero_and_preserves_existing_log(self):
        run = self.root / "make"
        (run / "logs").mkdir(parents=True)
        log = run / "logs/check_cxx.log"
        log.write_text("existing scientific check log\n")
        failing = run / "fail.py"
        failing.write_text("import sys\nprint('new failed run')\nsys.exit(7)\n")
        makefile = pathlib.Path(check_cxx.__file__).resolve().parents[1] / "Makefile"
        result = subprocess.run(["make", "--no-print-directory", "-f", str(makefile), "cxx",
                                 f"PY={sys.executable} {failing}"], cwd=run, text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(log.read_text(), "existing scientific check log\n")
        self.assertIn("new failed run", result.stderr)


if __name__ == "__main__":
    unittest.main()
