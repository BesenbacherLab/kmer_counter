"""Development tasks."""

import os
import re
import sys
from pathlib import Path
from shutil import which
from typing import List, Optional, Pattern

from duty import duty
from git_changelog.build import Changelog, Version
from jinja2.sandbox import SandboxedEnvironment

PY_SRC_PATHS = (Path(_) for _ in ("src", "tests", "duties.py"))
PY_SRC_LIST = tuple(str(_) for _ in PY_SRC_PATHS)
PY_SRC = " ".join(PY_SRC_LIST)
TESTING = os.environ.get("TESTING", "0") in {"1", "true"}
CI = os.environ.get("CI", "0") in {"1", "true", "yes", ""}
WINDOWS = os.name == "nt"
PTY = not WINDOWS and not CI


def latest(lines: List[str], regex: Pattern) -> Optional[str]:
    """
    Return the last released version.

    Arguments:
        lines: Lines of the changelog file.
        regex: A compiled regex to find version numbers.

    Returns:
        The last version.
    """
    for line in lines:
        match = regex.search(line)
        if match:
            return match.groupdict()["version"]
    return None


def unreleased(versions: List[Version], last_release: str) -> List[Version]:
    """
    Return the most recent versions down to latest release.

    Arguments:
        versions: All the versions (released and unreleased).
        last_release: The latest release.

    Returns:
        A list of versions.
    """
    for index, version in enumerate(versions):
        if version.tag == last_release:
            return versions[:index]
    return versions


@duty(silent=True)
def clean(ctx):
    """
    Delete temporary files.

    Arguments:
        ctx: The context instance (passed automatically).
    """
    ctx.run("rm -rf .coverage*")
    ctx.run("rm -rf .pytest_cache")
    ctx.run("rm -rf tests/.pytest_cache")
    ctx.run("rm -rf build")
    ctx.run("rm -rf dist")
    ctx.run("rm -rf pip-wheel-metadata")
    ctx.run("find . -type d -name __pycache__ | xargs rm -rf")
    ctx.run("find . -name '*.rej' -delete")


@duty
def release(ctx, version):
    """
    Release a new Python package.

    Arguments:
        ctx: The context instance (passed automatically).
        version: The new version number to use.
    """
    ctx.run(f"poetry version {version}", title=f"Bumping version in pyproject.toml to {version}", pty=PTY)
    ctx.run(f'sed -i "" "s/__version__ = .*/__version__ = \\"{version}\\"/g" src/kmer_counter/__init__.py')
    ctx.run("git add pyproject.toml src/kmer_counter/__init__.py", title="Staging files", pty=PTY)
    ctx.run(["git", "commit", "-m", f"chore: Prepare release {version}"], title="Committing changes", pty=PTY)
    ctx.run(f"git tag {version}", title="Tagging commit", pty=PTY)
    if not TESTING:
        ctx.run("git push", title="Pushing commits", pty=False)
        ctx.run("git push --tags", title="Pushing tags", pty=False)
        ctx.run("poetry build", title="Building dist/wheel", pty=PTY)
        ctx.run("poetry publish", title="Publishing version", pty=PTY)


@duty(silent=True)
def coverage(ctx):
    """
    Report coverage as text and HTML.

    Arguments:
        ctx: The context instance (passed automatically).
    """
    ctx.run("coverage combine .coverage-*", nofail=True)
    ctx.run("coverage report --rcfile=config/coverage.ini", capture=False)
    ctx.run("coverage html --rcfile=config/coverage.ini")


@duty
def test(ctx, match: str = ""):
    """
    Run the test suite.

    Arguments:
        ctx: The context instance (passed automatically).
        match: A pytest expression to filter selected tests.
    """
    py_version = f"{sys.version_info.major}{sys.version_info.minor}"
    os.environ["COVERAGE_FILE"] = f".coverage-{py_version}"
    ctx.run(
        ["pytest", "-c", "config/pytest.ini", "-n", "auto", "-k", match, "tests"],
        title="Running tests",
        pty=PTY,
    )
