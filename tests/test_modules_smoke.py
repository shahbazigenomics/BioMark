"""Smoke tests: every benchmark module without dedicated tests runs end to end
and returns the structure that main.py and report.py rely on.

These check structure and robustness, not score values, so they pass on any hardware.
"""
import importlib
import warnings
from unittest.mock import patch

import pytest

MODULES = [
    ("biomark.modules.RNA.scrna_benchmark", "run_scrna_benchmark"),
    ("biomark.modules.Protein.protein_benchmark", "run_protein_benchmark"),
    ("biomark.modules.Epigenomics.epigenomics_benchmark", "run_epigenomics_benchmark"),
    ("biomark.modules.Metagenomics.metagenomics_benchmark", "run_metagenomics_benchmark"),
    ("biomark.modules.LongRead.longread_benchmark", "run_longread_benchmark"),
]
IDS = [func for _, func in MODULES]

REQUIRED_KEYS = {"module", "score", "capability", "pipeline_steps", "total_time_seconds"}

# Deliberately weak machine: the benchmarks should report "cannot run", not crash.
TINY_MACHINE = {
    "total_ram_gb": 2.0,
    "available_ram_gb": 1.0,
    "total_ssd_gb": 20.0,
    "free_ssd_gb": 5.0,
    "cpu_cores": 2,
}

# Deliberately strong machine: forces every pipeline step to run, whatever host runs the tests.
BIG_MACHINE = {
    "total_ram_gb": 64.0,
    "available_ram_gb": 60.0,
    "total_ssd_gb": 2000.0,
    "free_ssd_gb": 1500.0,
    "cpu_cores": 16,
}


def _check_result(result):
    assert isinstance(result, dict)
    assert REQUIRED_KEYS <= set(result), f"missing keys: {REQUIRED_KEYS - set(result)}"
    assert isinstance(result["score"], (int, float))
    assert 0 <= result["score"] <= 100
    assert isinstance(result["capability"], str) and result["capability"]
    assert isinstance(result["pipeline_steps"], (dict, list))
    assert result["total_time_seconds"] >= 0


@pytest.mark.parametrize("module_name, func_name", MODULES, ids=IDS)
def test_module_runs_on_this_machine(module_name, func_name):
    module = importlib.import_module(module_name)
    _check_result(getattr(module, func_name)())


@pytest.mark.parametrize("module_name, func_name", MODULES, ids=IDS)
def test_module_survives_tiny_machine(module_name, func_name):
    module = importlib.import_module(module_name)
    with patch.object(module, "get_system_resources", return_value=dict(TINY_MACHINE)):
        _check_result(getattr(module, func_name)())


@pytest.mark.parametrize("module_name, func_name", MODULES, ids=IDS)
def test_module_survives_frozen_clock(module_name, func_name):
    """A step that finishes instantly must not cause a division by zero."""
    module = importlib.import_module(module_name)
    with patch("time.time", return_value=1000.0):
        _check_result(getattr(module, func_name)())


@pytest.mark.parametrize("module_name, func_name", MODULES, ids=IDS)
def test_module_runs_full_pipeline_on_big_machine(module_name, func_name):
    """On ample hardware the steps really run, so the score must be above zero."""
    module = importlib.import_module(module_name)
    with patch.object(module, "get_system_resources", return_value=dict(BIG_MACHINE)):
        result = getattr(module, func_name)()
    _check_result(result)
    assert result["score"] > 0


@pytest.mark.parametrize("module_name, func_name", MODULES, ids=IDS)
def test_big_machine_with_frozen_clock(module_name, func_name):
    """Every step executed with an instant clock must not divide by zero."""
    module = importlib.import_module(module_name)
    with patch.object(module, "get_system_resources", return_value=dict(BIG_MACHINE)), \
         patch("time.time", return_value=1000.0):
        _check_result(getattr(module, func_name)())


@pytest.mark.parametrize("module_name, func_name", MODULES, ids=IDS)
def test_no_numeric_warnings_leak_to_users(module_name, func_name):
    """Clients should never see NumPy RuntimeWarnings (e.g. spurious matmul warnings on macOS)."""
    module = importlib.import_module(module_name)
    with patch.object(module, "get_system_resources", return_value=dict(BIG_MACHINE)), \
         warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        _check_result(getattr(module, func_name)())
