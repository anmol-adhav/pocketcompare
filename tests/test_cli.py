"""CLI tests using pytest and subprocess."""
import subprocess
import sys
import pytest


def run(*args):
    return subprocess.run(
        [sys.executable, "-m", "pocketcompare.cli", *args],
        capture_output=True, text=True
    )


def test_cli_no_args():
    r = run()
    assert r.returncode != 0


def test_cli_help():
    r = run("--help")
    assert r.returncode == 0
    assert "homolog" in r.stdout or "usage" in r.stdout.lower()


def test_cli_homolog_help():
    r = run("homolog", "--help")
    assert r.returncode == 0
    assert "--pdb1" in r.stdout


def test_cli_mutant_help():
    r = run("mutant", "--help")
    assert r.returncode == 0
    assert "--wt" in r.stdout


@pytest.mark.integration
def test_cli_homolog_run():
    r = run(
        "homolog", "--pdb1", "2ZCP", "--pdb2", "5IYS",
        "--pairs", "A:172,A:173;A:176,A:177;A:225,A:237"
    )
    assert r.returncode == 0
    assert "RECOMMENDATION" in r.stdout
    assert "RMSD" in r.stdout


@pytest.mark.integration
def test_cli_mutant_self():
    r = run(
        "mutant", "--wt", "2ZCP", "--mut", "2ZCP",
        "--pocket", "A:172;A:176;A:225;A:228"
    )
    assert r.returncode == 0
    assert "UNCHANGED" in r.stdout
