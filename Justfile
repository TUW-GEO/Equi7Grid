# Default command lists all available recipes
_default:
    @just --list --unsorted

alias b := bump
alias c := clean
alias d := dist
alias h := hooks
alias q := check
alias t := test

PROJ := `uv version --short`

# lint python code using ruff
[private]
check-lint:
    uv run ruff check . --fix

# format python code using ruff
[private]
check-format:
    uv run ruff format .

# run the type checker ty
[private]
check-types:
    uv run ty check

# lint, format with ruff and type-check with ty
check: check-format check-lint check-types

# run tests with coverage
test:
    uv run pytest tests/

# run tests for all the supported Python versions
testall:
    uv run --python=3.11 pytest
    uv run --python=3.12 pytest
    uv run --python=3.13 pytest
    uv run --python=3.14 pytest


# run all the formatting, linting, and testing commands
ci PYTHON="3.12":
    uv run --python={{ PYTHON }} ruff format .
    uv run --python={{ PYTHON }} ruff check . --fix
    uv run --python={{ PYTHON }} ty check .
    uv run --python={{ PYTHON }} pytest tests/

# setup the pre-commit hooks
hooks:
    uvx pre-commit install

# print the current status of the project
status:
    @echo "Project Version: {{ PROJ }}"
    @echo "Running on: `uname`"

# clean all python build/compilation files and directories
clean: clean-build clean-pyc clean-test

# remove build artifacts
[private]
clean-build:
    rm -fr build/
    rm -fr _build/
    rm -fr dist/
    rm -fr .eggs/
    find . -name '*.egg-info' -exec rm -fr {} +
    find . -name '*.egg' -exec rm -f {} +

# remove Python file artifacts
[private]
clean-pyc:
    find . -name '*.pyc' -exec rm -f {} +
    find . -name '*.pyo' -exec rm -f {} +
    find . -name '*~' -exec rm -f {} +
    find . -name '__pycache__' -exec rm -fr {} +

# remove test and coverage artifacts
[private]
clean-test:
    rm -f .coverage
    rm -fr htmlcov/
    rm -fr .pytest_cache

# install dependencies in local venv
venv:
    uv sync

#[confirm("Do you really want to bump? (y/n)")]
[private]
prompt-confirm:

# bump the version, commit and add a tag <major|minor|patch|...>
bump INCREMENT="patch": && tag
    @uv version --bump {{ INCREMENT }} --dry-run
    @just prompt-confirm
    uv version --bump {{ INCREMENT }}

# tag the latest version
tag VERSION=`uv version --short`:
    git add pyproject.toml
    git add uv.lock
    git commit -m "Bumped version to {{VERSION}}"
    git tag -a "v{{VERSION}}"

# preview the documentation locally (serve the myst website)
docs:
    uv run myst

# build the source distribution and wheel file
dist:
    uv build
