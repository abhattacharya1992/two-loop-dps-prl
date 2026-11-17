"""Tests of the anc/ selection of make_arxiv.py (no TeX run).

The selection rule is tested on synthetic path lists; the archive layout is tested on this
repository when it is a git clone (skipped otherwise, e.g. in the extracted ancillary files).
"""
import shutil
import subprocess
import tarfile
import tempfile
import unittest
from pathlib import Path

import make_arxiv as ma

PREFIX = "paper2/manuscript_prl"


def in_git_clone() -> bool:
    if shutil.which("git") is None:
        return False
    run = subprocess.run(["git", "-C", str(ma.HERE), "rev-parse", "--is-inside-work-tree"],
                         capture_output=True, text=True)
    return run.returncode == 0 and run.stdout.strip() == "true"


class Selection(unittest.TestCase):
    TRACKED = [
        "Makefile", "README.md", "common/dps/results_hh.json", "paper2/hh/hh_two_loop/points/a.json",
        f"{PREFIX}/letter.tex", f"{PREFIX}/supplement.tex", f"{PREFIX}/tables/table_r0.tex",
        f"{PREFIX}/figures/dps_diagram.pdf", f"{PREFIX}/figures/dps_diagram.tex",
        f"{PREFIX}/figures/dps_diagram_tikz.tex",
        f"{PREFIX}/letter.pdf", f"{PREFIX}/supplement.pdf", f"{PREFIX}/make_tables.py",
        f"{PREFIX}/release/PRL_DPS_arxiv_sources.tar.gz",
        "derivations/README.md", "derivations/detailed_derivations.tex",
        "derivations/detailed_derivations.pdf", "derivations/checks/extension_results.json", "",
    ]
    TOP = {"letter.tex", "supplement.tex", "tables/table_r0.tex", "figures/dps_diagram.pdf"}

    def test_arxiv_only_supplemental_material_citation(self):
        original = ("See Supplemental Material at [URL will be inserted by publisher] "
                    "for the two-loop DPS derivation.")
        arxiv = ("See Supplemental Material accompanying this arXiv submission "
                 "for the two-loop DPS derivation.")
        for name in ("sm.bib", "letter.bbl"):
            if name == "letter.bbl":
                raw = original.replace("be inserted", "be\n  inserted").encode()
            else:
                raw = original.encode()
            self.assertEqual(ma.arxiv_member_bytes(Path(name), raw).decode(), arxiv)
        self.assertEqual(ma.arxiv_member_bytes(Path("supplement.tex"), original.encode()),
                         original.encode())
        with self.assertRaisesRegex(RuntimeError, "Expected one APS"):
            ma.arxiv_member_bytes(Path("sm.bib"), b"Reference already edited")

    def test_top_level_files_pdfs_and_release_are_left_out(self):
        paths, dirty = ma.select_ancillary(self.TRACKED, [], PREFIX, self.TOP)
        self.assertEqual(dirty, [])
        self.assertEqual(paths, sorted([
            "Makefile", "README.md", "common/dps/results_hh.json", "paper2/hh/hh_two_loop/points/a.json",
            f"{PREFIX}/make_tables.py"]))

    def test_unsafe_helper_and_old_model_copies_are_excluded(self):
        names = ['compile_manuscipt.sh', 'manuscipt_compile.sh',
                 'paper2/kernel/model_exact2_before_z6p_fix.py',
                 'paper2/kernel/model_fast_before_z6p_fix.py',
                 'paper2/kernel/model_exact2.py', 'paper2/kernel/model_fast.py']
        selected, dirty = ma.select_ancillary(names, names, PREFIX, set())
        self.assertEqual(selected, ['paper2/kernel/model_exact2.py',
                                    'paper2/kernel/model_fast.py'])
        self.assertEqual(dirty, selected)

    def test_repository_notes_are_excluded_including_changed_files(self):
        tracked = ["derivations/README.md", "derivations/detailed_derivations.tex",
                   "derivations/detailed_derivations.pdf", "derivations/checks/validation_runs.json",
                   "derivations-extra/keep.txt", "checks/README.md"]
        selected, dirty = ma.select_ancillary(tracked, tracked, PREFIX, self.TOP)
        self.assertEqual(selected, ["checks/README.md", "derivations-extra/keep.txt"])
        self.assertEqual(dirty, selected)

    def test_reference_diagrams_are_excluded_from_arxiv_archive(self):
        reference_files = [
            f"{PREFIX}/figures/{name}"
            for name in ("README.md", "collinear_diagrams_qq.tex", "collinear_diagrams_qq.pdf",
                         "collinear_diagrams_gg.tex", "collinear_diagrams_gg.pdf",
                         "selfenergy_soft_diagrams.tex", "selfenergy_soft_diagrams.pdf",
                         "r0_diagrams_captions.tex")
        ]
        selected, changed = ma.select_ancillary(reference_files, reference_files, PREFIX, self.TOP)
        self.assertEqual(selected, [])
        self.assertEqual(changed, [])

    def test_changes_count_only_for_files_that_go_into_anc(self):
        changed = [f"{PREFIX}/letter.pdf", f"{PREFIX}/tables/table_r0.tex",
                   f"{PREFIX}/figures/dps_diagram.tex",
                   f"{PREFIX}/release/PRL_DPS_arxiv_sources.tar.gz", "common/dps/results_hh.json",
                   "derivations/README.md", "derivations/detailed_derivations.tex"]
        _, dirty = ma.select_ancillary(self.TRACKED, changed, PREFIX, self.TOP)
        self.assertEqual(dirty, ["common/dps/results_hh.json"])

    def test_tex_sources_are_excluded_from_ancillary_files(self):
        tracked = [f"{PREFIX}/figures/dps_diagram.tex",
                   f"{PREFIX}/figures/dps_diagram_tikz.tex",
                   "paper2/hh/letter.tex", "paper2/hh/README.md"]
        paths, dirty = ma.select_ancillary(tracked, tracked, PREFIX, set())
        self.assertEqual(paths, ["paper2/hh/README.md"])
        self.assertEqual(dirty, paths)


@unittest.skipUnless(in_git_clone(), "not a git clone")
class Archive(unittest.TestCase):
    def test_layout(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp) / "arxiv.tar.gz"
            _, top_count, anc_count = ma.build_archive(dest, with_anc=True, allow_dirty=True)
            ma.check_layout(dest, with_anc=True)
            with tarfile.open(dest, "r:gz") as bundle:
                names = bundle.getnames()
                readme = bundle.extractfile(f"anc/{ma.ANC_README}").read().decode()
                arxiv_citations = {
                    member: bundle.extractfile(member).read().decode("utf-8")
                    for member in ("sm.bib", "letter.bbl")
                }
        top = [n for n in names if not n.startswith("anc/")]
        anc = [n for n in names if n.startswith("anc/")]
        self.assertEqual((len(top), len(anc)), (top_count, anc_count))
        self.assertIn("letter.tex", top)
        self.assertIn("supplement.tex", top)
        for note in arxiv_citations.values():
            self.assertIn("accompanying this arXiv submission", note)
            self.assertNotIn("[URL will be", note)
        self.assertIn("anc/common/dps/results_hh.json", anc)
        self.assertIn(f"anc/{PREFIX}/make_tables.py", anc)
        self.assertFalse(any(n.startswith("anc/derivations/") or n.startswith("derivations/") for n in names))
        self.assertNotIn(f"anc/{PREFIX}/letter.tex", anc)
        self.assertFalse(any(n.endswith(".tex") for n in anc))
        self.assertFalse(any(n.startswith(f"anc/{PREFIX}/release/") for n in anc))
        self.assertIn("at commit", readme)

    def test_archive_layout_rejects_repository_notes(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp) / "arxiv.tar.gz"
            with tarfile.open(dest, "w:gz") as bundle:
                for name in ("letter.tex", "anc/README_ANCILLARY.txt", "anc/derivations/README.md"):
                    payload = b"test"
                    member = tarfile.TarInfo(name)
                    member.size = len(payload)
                    bundle.addfile(member, __import__('io').BytesIO(payload))
            with self.assertRaisesRegex(RuntimeError, 'derivations'):
                ma.check_layout(dest, with_anc=True)

    def test_archive_layout_rejects_tex_in_ancillary_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp) / "arxiv.tar.gz"
            with tarfile.open(dest, "w:gz") as bundle:
                for name in ("letter.tex", "anc/README_ANCILLARY.txt",
                             "anc/paper2/manuscript_prl/figures/dps_diagram.tex"):
                    payload = b"test"
                    member = tarfile.TarInfo(name)
                    member.size = len(payload)
                    bundle.addfile(member, __import__('io').BytesIO(payload))
            with self.assertRaisesRegex(RuntimeError, "TeX files"):
                ma.check_layout(dest, with_anc=True)

    def test_without_anc(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp) / "arxiv.tar.gz"
            _, _, anc_count = ma.build_archive(dest, with_anc=False)
            ma.check_layout(dest, with_anc=False)
            with tarfile.open(dest, "r:gz") as bundle:
                self.assertFalse(any(n.startswith("anc/") for n in bundle.getnames()))
        self.assertEqual(anc_count, 0)


if __name__ == "__main__":
    unittest.main()
