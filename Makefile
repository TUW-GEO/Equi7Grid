.ONESHELL:
SHELL = /bin/bash

.PHONY: help clean environment install test version dist

CONDA_ENV_DIR = $(shell conda info --base)/envs/equi7grid
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate
TEST_ENV =

help:
	@echo "make clean"
	@echo " clean all python build/compilation files and directories"
	@echo "make venv"
	@echo " create the base virtualenv environment for the project"
	@echo "make conda"
	@echo " create the base conda environment for the project"
	@echo "make test"
	@echo " run test with coverage"
	@echo "make version"
	@echo " update _version.py with current version tag"
	@echo "make dist"
	@echo " build the package ready for distribution and update the version tag"

clean:
	find . -name '*.pyc' -exec rm --force {} +
	find . -name '*.pyo' -exec rm --force {} +
	find . -name '*~' -exec rm --force {} +
	rm --force .coverage
	rm --force --recursive .pytest_cache
	rm --force --recursive build/
	rm --force --recursive dist/
	rm --force --recursive *.egg-info

$(CONDA_ENV_DIR):
	@echo "creating new base equi7grid conda environment..."
	conda create -y -c conda-forge -n equi7grid python=3.10 pip mamba
	$(CONDA_ACTIVATE) equi7grid
	mamba install -y -c conda-forge gdal numpy scipy
	if [ $(TEST_ENV) ]; then pip install -e .[testing]; else pip install -e .; fi
	@echo "... finished."

conda: $(CONDA_ENV_DIR)
	@echo -e "conda environment is ready. To activate use:\n\tconda activate equi7grid"

venv/bin/activate:
	@echo "creating new base equi7grid virtualenv environment..."
	python3 -m venv venv
	source venv/bin/activate; pip install pygdal=="$(shell gdal-config --version).*"; if [ $(TEST_ENV) ]; then pip install -e .[testing]; else pip install -e .; fi
	@echo "... finished."

venv: venv/bin/activate
	@echo -e "virtualenv environment is ready. To activate use:\n\tsource venv/bin/activate"

test:
	pytest -rsx --verbose --color=yes --cov=equi7grid --cov-report term-missing

version:
	@version=$(shell git describe --always --tags | awk -F'[-]' '{if (NF == 1) print $$1; else printf "%s.dev%s+g%s\n", $$1, $$2, substr($$3, 2)}')
	echo -e "__version__ = \"$$version\"\n__commit__ = \"$(shell git rev-parse --short HEAD)\"" > src/equi7grid/_version.py

dist: version
	pip3 install build twine
	python3 -m build

