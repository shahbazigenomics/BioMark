# BioMark - RNA-seq Module Tests
# Author: Amir Shahbazi

import pytest
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from modules.RNA.rnaseq_benchmark import (
    generate_rna_sequence,
    generate_quality_scores,
    get_system_resources,
    evaluate_hardware_rnaseq,
    simulate_fastqc_multiqc,
    simulate_trimming,
    simulate_post_processing,
    simulate_featurecounts,
    simulate_salmon,
)


def test_generate_rna_sequence():
    seq = generate_rna_sequence(100)
    assert len(seq) == 100
    assert all(b in "ATGC" for b in seq)


def test_get_system_resources():
    res = get_system_resources()
    assert "total_ram_gb" in res
    assert res["total_ram_gb"] > 0


def test_evaluate_hardware_rnaseq():
    resources = get_system_resources()
    warnings  = evaluate_hardware_rnaseq(resources)
    assert isinstance(warnings, list)
    assert len(warnings) > 0


def test_simulate_fastqc_multiqc():
    result = simulate_fastqc_multiqc(n_reads=600, n_samples=6)
    assert "samples_analyzed" in result
    assert result["samples_analyzed"] == 6


def test_simulate_trimming():
    result = simulate_trimming(n_reads=600, n_samples=6)
    assert "tool_comparison" in result
    assert "Trimmomatic" in result["tool_comparison"]
    assert "fastp" in result["tool_comparison"]


def test_simulate_post_processing():
    result = simulate_post_processing(n_reads=1000, n_samples=6)
    assert "status" in result
    assert "substeps" in result


def test_simulate_featurecounts():
    result = simulate_featurecounts(n_reads=1000, n_samples=6)
    assert "genes_in_annotation" in result
    assert result["samples_quantified"] == 6


def test_simulate_salmon():
    result = simulate_salmon(n_reads=1000, n_samples=6)
    assert "transcripts_quantified" in result
    assert len(result["sample_details"]) == 6
