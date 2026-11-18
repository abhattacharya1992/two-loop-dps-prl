"""Regression checks for the frozen-commit comparison harness."""
from contextlib import contextmanager
import io
import json
from pathlib import Path
import subprocess
import tarfile
import tempfile
import unittest
from unittest.mock import patch

import release_compare as rc


@contextmanager
def directory(path):
    before = Path.cwd()
    import os
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(before)


def commit(root, content):
    (root / "x.txt").write_text(content)
    subprocess.run(["git", "add", "x.txt"], cwd=root, check=True, capture_output=True)
    subprocess.run(["git", "commit", "-m", "test"], cwd=root, check=True, capture_output=True)
    return rc.git(root, "rev-parse", "HEAD")


class ComparisonTests(unittest.TestCase):
    def test_science_lines_normalize_only_temp_paths(self):
        s = """PASS: 66 checks at /tmp/tmp.Ab09/stage1.json
r0_theory_quarks.json reproduced, largest difference 0.0e+00
A regular output that should not be compared
"""
        self.assertEqual(rc.science_lines(s), [
            "PASS: 66 checks at /tmp/<temporary>/stage1.json",
            "r0_theory_quarks.json reproduced, largest difference 0.0e+00",
        ])
        changed = s.replace("66 checks", "65 checks")
        self.assertNotEqual(rc.science_lines(s), rc.science_lines(changed))

    def test_blob_inventory_reads_committed_file_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            subprocess.run(["git", "init", "-q"], cwd=root, check=True)
            subprocess.run(["git", "config", "user.name", "Researcher"], cwd=root, check=True)
            subprocess.run(["git", "config", "user.email", "test@example.org"], cwd=root, check=True)
            revision = commit(root, "content")
            tree = rc.git_blobs(root, revision)
            self.assertEqual(len(tree), 1)
            self.assertEqual(len(tree["x.txt"]), 40)

    def test_historical_cache_is_not_a_required_scientific_input(self):
        reference = {
            "paper2/r0/theory/stored/history/2026-10-22-pr1-initial-fdh/fdh_laurent_x1o2.json": "old",
            "paper2/r0/theory/r0_theory_quarks.json": "same",
            "paper2/manuscript_prl/letter.tex": "paper",
            "common/dps/dps/analysis.py": "previous",
        }
        candidate = {
            "paper2/r0/theory/r0_theory_quarks.json": "same",
            "paper2/manuscript_prl/letter.tex": "paper",
            "common/dps/dps/analysis.py": "revised",
        }
        checked, differences = rc.compare_inputs(reference, candidate)
        self.assertEqual(len(checked), 2)
        self.assertEqual(differences, [])
        candidate["paper2/r0/theory/r0_theory_quarks.json"] = "changed"
        self.assertEqual(len(rc.compare_inputs(reference, candidate)[1]), 1)
        del candidate["paper2/manuscript_prl/letter.tex"]
        self.assertEqual(len(rc.compare_inputs(reference, candidate)[1]), 2)

    def test_archive_compares_only_publication_top_level(self):
        with tempfile.TemporaryDirectory() as tmp:
            archive = Path(tmp) / "sources.tar.gz"
            with tarfile.open(archive, "w:gz") as bundle:
                for name, data in (("letter.tex", b"publication"),
                                   ("anc/README_ANCILLARY.txt", b"commit-specific"),
                                   ("anc/data.json", b"stored input")):
                    entry = tarfile.TarInfo(name)
                    entry.size = len(data)
                    bundle.addfile(entry, io.BytesIO(data))
            self.assertEqual(rc.archive_sources(archive), {"letter.tex": rc.sha256(b"publication")})

    def test_generated_files_are_compared_by_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("main", "clean"):
                (root / name / "paper2/manuscript_prl/tables").mkdir(parents=True)
                (root / name / "paper2/manuscript_prl/numbers.tex").write_text("123")
                (root / name / "paper2/manuscript_prl/tables/table_zz.tex").write_text("abcd")
            a, b = root / "main", root / "clean"
            paths = [Path("paper2/manuscript_prl/numbers.tex"),
                     Path("paper2/manuscript_prl/tables/table_zz.tex")]
            self.assertEqual(rc.compare_files(a, b, paths), [])
            (b / paths[1]).write_text("abce")
            self.assertEqual(rc.compare_files(a, b, paths), [str(paths[1])])

    def test_git_worktrees_are_temporary_and_leave_branch_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            subprocess.run(["git", "init", "-q"], cwd=root, check=True)
            subprocess.run(["git", "config", "user.name", "Researcher"], cwd=root, check=True)
            subprocess.run(["git", "config", "user.email", "test@example.org"], cwd=root, check=True)
            baseline = commit(root, "original")
            candidate = commit(root, "revised")
            expected_branch = rc.git(root, "symbolic-ref", "--short", "HEAD")
            with rc.worktrees(root, baseline, candidate) as (a, b):
                self.assertEqual((a / "x.txt").read_text(), "original")
                self.assertEqual((b / "x.txt").read_text(), "revised")
            self.assertEqual(rc.git(root, "rev-parse", "HEAD"), candidate)
            self.assertEqual(rc.git(root, "symbolic-ref", "--short", "HEAD"), expected_branch)
            self.assertEqual(len(rc.git(root, "worktree", "list").splitlines()), 1)

    def test_quick_mode_creates_a_report_without_mutating_repository(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            subprocess.run(["git", "init", "-q"], cwd=root, check=True)
            subprocess.run(["git", "config", "user.name", "Researcher"], cwd=root, check=True)
            subprocess.run(["git", "config", "user.email", "test@example.org"], cwd=root, check=True)
            baseline = commit(root, "one")
            candidate = commit(root, "two")
            with directory(root):
                self.assertEqual(rc.main(["--baseline", baseline, "--candidate", candidate,
                                          "--output", "out", "--quick"]), 0)
            report = json.loads((root / "out/report.json").read_text())
            self.assertTrue(report["complete"])
            self.assertEqual(report["mode"], "quick")
            self.assertIn("calculations not run", (root / "out/report.md").read_text())
            self.assertEqual(rc.git(root, "rev-parse", "HEAD"), candidate)

    def test_result_is_not_certified_on_failed_inputs(self):
        result = {"baseline": "a", "candidate": "b", "mode": "full",
                  "input_files_checked": 5, "input_mismatches": [1],
                  "failures": ["one scientific input differs"]}
        self.assertIn("FAILED", rc.markdown(result))
        self.assertIn("one scientific input differs", rc.markdown(result))


if __name__ == "__main__":
    unittest.main()
