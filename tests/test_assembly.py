# BioMark - Assembly Module Tests
# Author: Amir Shahbazi

import pytest
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from modules.Assembly.assembly_benchmark import (
    get_system_resources,
    evaluate_hardware,
    simulate_spades,
    simulate_flye,
    simulate_hifiasm,
    simulate_quast,
    simulate_busco,
    run_assembly_benchmark,
)


def test_get_system_resources():
    res = get_system_resources()
    assert "total_ram_gb" in res
    assert res["total_ram_gb"] > 0


def test_evaluate_hardware():
    resources = get_system_resources()
    warnings  = evaluate_hardware(resources)
    assert isinstance(warnings, list)
    assert len(warnings) > 0


def test_simulate_spades():
    resources = get_system_resources()
    result    = simulate_spades(resources)
    assert "step" in result
    assert "status" in result
    assert result["step"] == "SPAdes"


def test_simulate_flye():
    resources = get_system_resources()
    result    = simulate_flye(resources)
    assert result["step"] == "Flye"


def test_simulate_hifiasm():
    resources = get_system_resources()
    result    = simulate_hifiasm(resources)
    assert result["step"] == "Hifiasm"


def test_simulate_quast():
    resources = get_system_resources()
    result    = simulate_quast(resources)
    assert "step" in result
    assert "status" in result


def test_simulate_busco():
    resources = get_system_resources()
    result    = simulate_busco(resources)
    assert "step" in result
    assert "status" in result


def test_run_assembly_benchmark():
    result = run_assembly_benchmark()
    assert result["module"] == "Assembly"
    assert 0 <= result["score"] <= 100
    assert "pipeline_steps" in result
