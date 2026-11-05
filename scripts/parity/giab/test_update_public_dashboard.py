#!/usr/bin/env python3
"""Focused tests for Pages history merge and workflow provenance.

No GitHub contact. Git fixtures are local bare repositories.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[3]
UPDATER = REPO / "scripts" / "parity" / "giab" / "update_public_dashboard.py"
PUBLISHER = REPO / "scripts" / "ci" / "publish_parity_site.sh"


def pages_history() -> dict:
    return {
        "meta": {
            "title": "keep-me",
            "notes": "original-notes",
            "extra": "stay",
            "updated_utc": "2026-08-24T09:40:53Z",
        },
        "runs": [
            {
                "id": "pages-only",
                "workflow": "nightly-equivalence",
                "commit_sha": "7611aec4df1ae21e300af62dc5e2ebc801421465",
                "metrics": [{"f1": None, "region": "chr20"}],
                "scope": {"regions": ["chr20"], "samples": ["HG002", "HG003", "HG004"]},
            },
            {
                "id": "genomewide-abc123-ci-subset",
                "workflow": "genomewide-validation",
                "commit_sha": "abc123",
                "metrics": [{"f1": 0.01, "engine": "rust", "variant_type": "SNP"}],
            },
        ],
    }


def write_run_dir(path: Path, mode: str = "ci-subset") -> None:
    path.mkdir(parents=True, exist_ok=True)
    (path / "SCOPE.txt").write_text(f"MODE {mode}\n", encoding="utf-8")
    (path / "intervals.txt").write_text("20:1-100\n", encoding="utf-8")
    row = {
        "sample": "HG001",
        "mode": mode,
        "equiv_results": {
            "rust_vs_truth": [
                {
                    "query_label": "rust",
                    "stratum": "ALL",
                    "snp": {"precision": 0.2, "recall": 0.3, "f1": 0.4},
                    "indel": {"precision": 0.5, "recall": 0.6, "f1": 0.7},
                }
            ],
            "java_vs_truth": [
                {
                    "query_label": "java",
                    "stratum": "ALL",
                    "snp": {"precision": 0.2, "recall": 0.3, "f1": 0.4},
                    "indel": {"precision": 0.5, "recall": 0.6, "f1": 0.7},
                }
            ],
        },
    }
    (path / "samples.jsonl").write_text(json.dumps(row) + "\n", encoding="utf-8")


def run_updater(tmp: Path, *extra: str) -> subprocess.CompletedProcess[str]:
    site = tmp / "site"
    run_dir = tmp / "run"
    write_run_dir(run_dir)
    return subprocess.run(
        [
            "python3",
            str(UPDATER),
            "--source",
            "genomewide",
            "--run-dir",
            str(run_dir),
            "--site-dir",
            str(site),
            "--commit-sha",
            "abc123",
            *extra,
        ],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=False,
    )


class UpdatePublicDashboardTest(unittest.TestCase):
    def test_pages_only_record_survives_and_new_record_is_appended(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            tmp = Path(raw)
            pages = tmp / "pages.json"
            original = pages_history()
            # Drop the id this run would generate so the new row is actually new.
            original["runs"] = [original["runs"][0]]
            pages.write_text(json.dumps(original), encoding="utf-8")
            checkout = {"meta": {"title": "stale"}, "runs": [{"id": "checkout-only", "workflow": "stale"}]}
            (tmp / "checkout.json").write_text(json.dumps(checkout), encoding="utf-8")
            before = pages.read_bytes()
            result = run_updater(
                tmp,
                "--workflow",
                "giab-genomewide",
                "--pages-history-status",
                "present",
                "--pages-history",
                str(pages),
                "--fallback-history",
                str(tmp / "checkout.json"),
                "--max-runs",
                "1",
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(pages.read_bytes(), before)
            written = json.loads((tmp / "site" / "data" / "history.json").read_text(encoding="utf-8"))
            self.assertEqual([row["id"] for row in written["runs"]], ["pages-only", "genomewide-abc123-ci-subset"])
            self.assertEqual(written["runs"][0], original["runs"][0])
            self.assertNotIn("checkout-only", [row["id"] for row in written["runs"]])
            self.assertEqual(written["meta"]["title"], "keep-me")
            self.assertEqual(written["meta"]["notes"], "original-notes")
            self.assertEqual(written["meta"]["extra"], "stay")
            self.assertNotEqual(written["meta"]["updated_utc"], "2026-08-24T09:40:53Z")
            self.assertEqual(written["runs"][1]["workflow"], "giab-genomewide")

    def test_same_record_id_is_not_duplicated_or_rewritten(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            tmp = Path(raw)
            pages = tmp / "pages.json"
            original = pages_history()
            pages.write_text(json.dumps(original, indent=2) + "\n", encoding="utf-8")
            before = pages.read_bytes()
            result = run_updater(
                tmp,
                "--workflow",
                "giab-genomewide",
                "--pages-history-status",
                "present",
                "--pages-history",
                str(pages),
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("unchanged", result.stdout)
            self.assertEqual(pages.read_bytes(), before)
            self.assertFalse((tmp / "site" / "data" / "history.json").exists())
            self.assertFalse((tmp / "site" / "data" / "latest.json").exists())

    def test_workflow_name_comes_from_the_caller(self) -> None:
        for workflow in ("giab-genomewide", "genomewide-validation"):
            with tempfile.TemporaryDirectory() as raw:
                tmp = Path(raw)
                pages = tmp / "pages.json"
                pages.write_text(json.dumps({"meta": {}, "runs": []}), encoding="utf-8")
                result = run_updater(
                    tmp,
                    "--workflow",
                    workflow,
                    "--pages-history-status",
                    "present",
                    "--pages-history",
                    str(pages),
                    "--commit-sha",
                    "def456",
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                written = json.loads((tmp / "site" / "data" / "history.json").read_text(encoding="utf-8"))
                self.assertEqual(written["runs"][0]["workflow"], workflow)
                self.assertEqual(written["runs"][0]["id"], "genomewide-def456-ci-subset")

    def test_missing_pages_history_uses_checkout_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            tmp = Path(raw)
            fallback = {"meta": {"title": "seed"}, "runs": [{"id": "seed-row", "workflow": "joint-cohort-scale"}]}
            (tmp / "fallback.json").write_text(json.dumps(fallback), encoding="utf-8")
            result = run_updater(
                tmp,
                "--workflow",
                "giab-genomewide",
                "--pages-history-status",
                "absent",
                "--fallback-history",
                str(tmp / "fallback.json"),
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("checkout history", result.stderr)
            written = json.loads((tmp / "site" / "data" / "history.json").read_text(encoding="utf-8"))
            self.assertEqual([row["id"] for row in written["runs"]], ["seed-row", "genomewide-abc123-ci-subset"])
            self.assertEqual(written["runs"][0], fallback["runs"][0])
            self.assertEqual(written["meta"]["title"], "seed")

    def test_unavailable_or_malformed_history_is_not_replaced(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            tmp = Path(raw)
            sentinel = tmp / "site" / "data"
            sentinel.mkdir(parents=True)
            (sentinel / "history.json").write_text("SENTINEL", encoding="utf-8")
            unavailable = run_updater(
                tmp,
                "--workflow",
                "giab-genomewide",
                "--pages-history-status",
                "unavailable",
            )
            self.assertEqual(unavailable.returncode, 2)
            self.assertEqual((sentinel / "history.json").read_text(encoding="utf-8"), "SENTINEL")
            self.assertFalse((sentinel / "latest.json").exists())

            pages = tmp / "broken.json"
            pages.write_text("{", encoding="utf-8")
            malformed = run_updater(
                tmp,
                "--workflow",
                "genomewide-validation",
                "--pages-history-status",
                "present",
                "--pages-history",
                str(pages),
            )
            self.assertEqual(malformed.returncode, 2)
            self.assertIn("malformed", malformed.stderr)
            self.assertEqual((sentinel / "history.json").read_text(encoding="utf-8"), "SENTINEL")


def git(repo: Path, *args: str) -> None:
    subprocess.run(["git", "-C", str(repo), *args], check=True, capture_output=True, text=True)


def git_bytes(repo: Path, *args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(repo), *args],
        check=True,
        capture_output=True,
        text=True,
    ).stdout


class PublishParitySiteTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.bare = self.root / "pages.git"
        subprocess.run(["git", "init", "--bare", "-b", "gh-pages", str(self.bare)], check=True, capture_output=True)
        self.checkout = self.root / "checkout"
        self._init_checkout()

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def _init_checkout(self) -> None:
        self.checkout.mkdir()
        script_dir = self.checkout / "scripts" / "parity" / "giab"
        script_dir.mkdir(parents=True)
        os.symlink(UPDATER, script_dir / "update_public_dashboard.py")
        (self.checkout / "docs" / "parity-site" / "data").mkdir(parents=True)
        (self.checkout / "README").write_text("checkout\n", encoding="utf-8")
        write_run_dir(self.checkout / "run")
        git(self.checkout, "init", "-b", "main")
        git(self.checkout, "config", "user.email", "test@example.com")
        git(self.checkout, "config", "user.name", "test")
        git(self.checkout, "add", "README")
        git(self.checkout, "commit", "-m", "init")

    def _seed_pages(self, history: dict, index: str = "<p>site</p>\n") -> None:
        seed = self.root / "seed"
        subprocess.run(["git", "init", "-b", "gh-pages", str(seed)], check=True, capture_output=True)
        git(seed, "config", "user.email", "test@example.com")
        git(seed, "config", "user.name", "test")
        (seed / "index.html").write_text(index, encoding="utf-8")
        (seed / "data").mkdir()
        (seed / "data" / "history.json").write_text(json.dumps(history, indent=2) + "\n", encoding="utf-8")
        git(seed, "add", "index.html", "data/history.json")
        git(seed, "commit", "-m", "pages history")
        git(seed, "remote", "add", "origin", str(self.bare))
        git(seed, "push", "origin", "HEAD:gh-pages")

    def _publish(self, workflow: str, sha: str, remote: str | None = None) -> subprocess.CompletedProcess[str]:
        env = os.environ.copy()
        env.pop("PARITY_SITE_PAGES_HISTORY_STATUS", None)
        env.pop("PARITY_SITE_PAGES_HISTORY_FILE", None)
        env.update(
            {
                "PARITY_SITE_REPO_ROOT": str(self.checkout),
                "PARITY_SITE_SOURCE": "genomewide",
                "PARITY_SITE_WORKFLOW": workflow,
                "PARITY_SITE_RUN_DIR": str(self.checkout / "run"),
                "PARITY_SITE_GIT_REMOTE": remote if remote is not None else str(self.bare),
                "PARITY_SITE_PUSH_PAGES": "1",
                "PARITY_SITE_OMIT_FROM_SITE_DEPLOY": "1",
                "PARITY_SITE_COMMIT": "0",
                "PARITY_SITE_PUSH": "0",
                "GITHUB_SHA": sha,
            }
        )
        return subprocess.run(
            ["bash", str(PUBLISHER)],
            cwd=self.checkout,
            env=env,
            text=True,
            capture_output=True,
            check=False,
        )

    def _pages_history(self) -> dict:
        text = git_bytes(self.bare, "show", "gh-pages:data/history.json")
        return json.loads(text)

    def test_publish_preserves_pages_rows_and_is_idempotent(self) -> None:
        original = pages_history()
        original["runs"] = [original["runs"][0]]
        kept = json.loads(json.dumps(original["runs"][0]))
        self._seed_pages(original)
        first = self._publish("giab-genomewide", "abc123")
        self.assertEqual(first.returncode, 0, first.stderr + first.stdout)
        written = self._pages_history()
        self.assertEqual([row["id"] for row in written["runs"]], ["pages-only", "genomewide-abc123-ci-subset"])
        self.assertEqual(written["runs"][0], kept)
        self.assertEqual(written["runs"][1]["workflow"], "giab-genomewide")
        self.assertEqual(git_bytes(self.bare, "show", "gh-pages:index.html"), "<p>site</p>\n")
        self.assertFalse((self.checkout / "docs" / "parity-site" / "data" / "history.json").exists())

        second = self._publish("giab-genomewide", "abc123")
        self.assertEqual(second.returncode, 0, second.stderr + second.stdout)
        self.assertIn("already present", second.stdout)
        again = self._pages_history()
        self.assertEqual(again["runs"], written["runs"])

        third = self._publish("genomewide-validation", "def456")
        self.assertEqual(third.returncode, 0, third.stderr + third.stdout)
        both = self._pages_history()
        self.assertEqual(
            [row["id"] for row in both["runs"]],
            ["pages-only", "genomewide-abc123-ci-subset", "genomewide-def456-ci-subset"],
        )
        self.assertEqual(both["runs"][0], kept)
        self.assertEqual(both["runs"][1], written["runs"][1])
        self.assertEqual(both["runs"][2]["workflow"], "genomewide-validation")
        self.assertEqual(git_bytes(self.bare, "show", "gh-pages:index.html"), "<p>site</p>\n")

    def test_absent_branch_uses_checkout_fallback(self) -> None:
        checkout_history = {
            "meta": {"title": "seed"},
            "runs": [{"id": "seed-row", "workflow": "joint-cohort-scale"}],
        }
        path = self.checkout / "docs" / "parity-site" / "data" / "history.json"
        path.write_text(json.dumps(checkout_history), encoding="utf-8")
        result = self._publish("giab-genomewide", "abc123")
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        self.assertIn("no gh-pages branch", result.stderr)
        written = self._pages_history()
        self.assertEqual([row["id"] for row in written["runs"]], ["seed-row", "genomewide-abc123-ci-subset"])
        self.assertEqual(written["runs"][0], checkout_history["runs"][0])

    def test_unavailable_remote_does_not_replace_checkout_history(self) -> None:
        path = self.checkout / "docs" / "parity-site" / "data" / "history.json"
        path.write_text('{"meta":{},"runs":[{"id":"local"}]}\n', encoding="utf-8")
        before = path.read_bytes()
        result = self._publish("giab-genomewide", "abc123", remote=str(self.root / "missing.git"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("refusing to publish", result.stderr)
        self.assertEqual(path.read_bytes(), before)

    def test_malformed_pages_history_is_left_in_place(self) -> None:
        seed = self.root / "seed-broken"
        subprocess.run(["git", "init", "-b", "gh-pages", str(seed)], check=True, capture_output=True)
        git(seed, "config", "user.email", "test@example.com")
        git(seed, "config", "user.name", "test")
        (seed / "data").mkdir()
        (seed / "data" / "history.json").write_text("{", encoding="utf-8")
        git(seed, "add", "data/history.json")
        git(seed, "commit", "-m", "broken")
        git(seed, "remote", "add", "origin", str(self.bare))
        git(seed, "push", "origin", "HEAD:gh-pages")
        result = self._publish("giab-genomewide", "abc123")
        self.assertEqual(result.returncode, 2, result.stderr + result.stdout)
        self.assertEqual(git_bytes(self.bare, "show", "gh-pages:data/history.json"), "{")


if __name__ == "__main__":
    unittest.main()
