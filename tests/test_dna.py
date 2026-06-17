# BioMark - DNA/WES Module Tests
# Author: Amir Shahbazi

import pytest
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from modules.DNA.wes_benchmark import (
    generate_dna_sequence,
    generate_quality_scores,
    simulate_fastqc,
    simulate_trimmomatic,
    simulate_variant_calling_gatk,
    simulate_variant_calling_bcftools,
    simulate_variant_calling_freebayes,
    simulate_annotation,
    simulate_filtering,
    get_system_resources,
    evaluate_hardware,
)


def test_generate_dna_sequence():
    seq = generate_dna_sequence(100)
    assert len(seq) == 100
    assert all(b in "ATGC" for b in seq)


def test_generate_quality_scores():
    qual = generate_quality_scores(100)
    assert len(qual) == 100
    assert all(33 <= ord(c) <= 74 for c in qual)


def test_get_system_resources():
    res = get_system_resources()
    assert "total_ram_gb" in res
    assert "free_ssd_gb" in res
    assert "cpu_cores" in res
    assert res["total_ram_gb"] > 0
    assert res["cpu_cores"] > 0


def test_evaluate_hardware():
    resources = get_system_resources()
    warnings  = evaluate_hardware(resources)
    assert isinstance(warnings, list)
    assert len(warnings) > 0
    for w in warnings:
        assert "level" in w
        assert "component" in w
        assert "message" in w


def test_simulate_fastqc():
    result = simulate_fastqc(n_reads=1000, read_length=150)
    assert "PASS" in result["status"]
    assert result["reads_analyzed"] == 1000
    assert "mean_gc" in result
    assert 0 <= result["mean_gc"] <= 100


def test_simulate_trimmomatic():
    result = simulate_trimmomatic(n_reads=1000, read_length=150)
    assert "PASS" in result["status"]
    assert result["reads_passed"] + result["reads_dropped"] == 1000
    assert 0 <= result["survival_rate"] <= 100


def test_simulate_variant_calling_gatk():
    result = simulate_variant_calling_gatk(n_positions=10000)
    assert "status" in result
    assert "variants_called" in result
    assert result["variants_called"] >= 0


def test_simulate_variant_calling_bcftools():
    result = simulate_variant_calling_bcftools(n_positions=10000)
    assert "status" in result
    assert "variants_called" in result


def test_simulate_variant_calling_freebayes():
    result = simulate_variant_calling_freebayes(n_positions=10000)
    assert "status" in result
    assert "variants_called" in result


def test_simulate_annotation():
    result = simulate_annotation(n_variants=1000)
    assert "status" in result
    assert "variants_annotated" in result
    assert result["variants_annotated"] == 1000


def test_simulate_filtering():
    result = simulate_filtering(n_variants=1000)
    assert "status" in result
    assert "variants_input" in result
    assert result["variants_input"] == 1000
    assert result["final_candidates"] <= 1000
