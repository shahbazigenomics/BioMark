# BioMark - Bulk RNA-seq Pipeline Benchmark Module
# Simulates real RNA-seq computational workloads
# Based on standard bulk RNA-seq pipeline
# Author: Amir Shahbazi
# GitHub: shahbazigenomics

import time
import random
import os
import gzip
import multiprocessing
import psutil
import numpy as np
from scipy import stats

# ════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ════════════════════════════════════════════════════════════

def generate_rna_sequence(length):
    """Generate RNA-like sequence (higher GC content than DNA)"""
    # RNA-seq reads tend to have higher GC content
    bases = ['A', 'T', 'G', 'C']
    weights = [0.22, 0.22, 0.28, 0.28]  # GC-biased
    return ''.join(random.choices(bases, weights=weights, k=length))

def generate_quality_scores(length):
    """Generate Phred quality scores"""
    return ''.join([chr(random.randint(33, 74)) for _ in range(length)])

def get_system_resources():
    mem  = psutil.virtual_memory()
    disk = psutil.disk_usage('/')
    return {
        "total_ram_gb"    : round(mem.total / (1024**3), 1),
        "available_ram_gb": round(mem.available / (1024**3), 1),
        "total_ssd_gb"    : round(disk.total / (1024**3), 1),
        "free_ssd_gb"     : round(disk.free / (1024**3), 1),
        "cpu_cores"       : multiprocessing.cpu_count(),
    }

def evaluate_hardware_rnaseq(resources):
    """
    Honest hardware evaluation for RNA-seq
    Different requirements from WES:
    - STAR needs ~32GB RAM for human genome index
    - HISAT2 needs ~8GB RAM
    - featureCounts/Salmon are less RAM intensive
    - DESeq2 in R needs ~8-16GB for large datasets
    """
    warnings = []
    ram   = resources["total_ram_gb"]
    ssd   = resources["free_ssd_gb"]
    cores = resources["cpu_cores"]

    # RAM Assessment for RNA-seq
    if ram < 8:
        warnings.append({
            "level"         : "🔴 CRITICAL",
            "component"     : "RAM",
            "message"       : f"{ram}GB RAM — Cannot run RNA-seq pipeline",
            "recommendation": "Minimum 8GB for HISAT2, 32GB for STAR alignment"
        })
    elif ram < 16:
        warnings.append({
            "level"         : "🔴 ERROR",
            "component"     : "RAM",
            "message"       : f"{ram}GB RAM — STAR alignment will fail (needs 32GB)",
            "recommendation": "Use HISAT2 instead of STAR, needs 8GB minimum"
        })
    elif ram < 32:
        warnings.append({
            "level"         : "🟡 WARNING",
            "component"     : "RAM",
            "message"       : f"{ram}GB RAM — HISAT2 ok, STAR may struggle",
            "recommendation": "32GB recommended for STAR + DESeq2 with large datasets"
        })
    else:
        warnings.append({
            "level"         : "🟢 GOOD",
            "component"     : "RAM",
            "message"       : f"{ram}GB RAM — Suitable for full RNA-seq pipeline",
            "recommendation": "Sufficient for STAR + DESeq2 + 3vs3 comparison"
        })

    # SSD Assessment for RNA-seq
    # FASTQ(10GB) + STAR index(32GB) + BAM(5-15GB) + results
    if ssd < 60:
        warnings.append({
            "level"         : "🔴 CRITICAL",
            "component"     : "SSD",
            "message"       : f"{ssd}GB free — Cannot store STAR index + samples",
            "recommendation": "Need 60GB+ free: STAR index(32GB) + FASTQ(10GB) + BAM(15GB)"
        })
    elif ssd < 200:
        warnings.append({
            "level"         : "🟡 WARNING",
            "component"     : "SSD",
            "message"       : f"{ssd}GB free — Very tight for 3vs3 RNA-seq",
            "recommendation": "200GB+ recommended for 6 samples (3 disease vs 3 control)"
        })
    else:
        warnings.append({
            "level"         : "🟢 GOOD",
            "component"     : "SSD",
            "message"       : f"{ssd}GB free — Sufficient for 3vs3 RNA-seq",
            "recommendation": "500GB+ for larger cohorts"
        })

    # CPU — only meaningful if RAM is sufficient for at least HISAT2
    ram_ok = ram >= 8
    if not ram_ok:
        warnings.append({
            "level"         : "⚪ N/A",
            "component"     : "CPU",
            "message"       : f"{cores} cores — CPU irrelevant while RAM is insufficient",
            "recommendation": "Fix RAM first — need 8GB+ minimum for HISAT2",
        })
    elif cores < 8:
        warnings.append({
            "level"         : "🟡 WARNING",
            "component"     : "CPU",
            "message"       : f"{cores} cores — STAR alignment will be slow",
            "recommendation": "8+ cores recommended, STAR uses all available cores"
        })
    else:
        warnings.append({
            "level"         : "🟢 GOOD",
            "component"     : "CPU",
            "message"       : f"{cores} cores — Good for RNA-seq pipeline",
            "recommendation": "16+ cores for faster STAR alignment"
        })

    return warnings


# ════════════════════════════════════════════════════════════
# STEP 1: QC SIMULATION (FastQC + MultiQC)
# ════════════════════════════════════════════════════════════

def simulate_fastqc_multiqc(n_reads=100000,
                             read_length=150,
                             n_samples=6):
    """
    Simulates FastQC + MultiQC workload for RNA-seq:
    - FastQC: per-sample QC (6 samples in 3vs3 design)
    - MultiQC: aggregate report across all samples
    - Check: per-base quality, GC content, adapter content,
             sequence duplication, overrepresented sequences
    Bottleneck: SSD read speed
    Real data: ~5-10GB FASTQ per sample × 6 samples
    """
    print("    📊 Step 1: QC Simulation (FastQC + MultiQC)...")
    start = time.time()

    sample_results = []

    for s in range(n_samples):
        sample_name = (
            f"Disease_{s+1}" if s < 3 else f"Control_{s-2}"
        )
        gc_distribution  = []
        qual_per_base    = [[] for _ in range(read_length)]
        dup_tracker      = {}
        adapter_count    = 0
        overrep_count    = 0

        for i in range(n_reads // n_samples):
            seq  = generate_rna_sequence(read_length)
            qual = generate_quality_scores(read_length)

            # GC content per read
            gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
            gc_distribution.append(gc)

            # Per-base quality
            for pos, q in enumerate(qual[:10]):
                qual_per_base[pos].append(ord(q) - 33)

            # Duplication check (simplified)
            key = seq[:20]
            dup_tracker[key] = dup_tracker.get(key, 0) + 1

            # Adapter detection (Illumina universal adapter)
            if 'AGATCGGAAGAGC' in seq:
                adapter_count += 1

            # Overrepresented sequences
            if dup_tracker.get(key, 0) > 5:
                overrep_count += 1

        mean_gc    = round(
            sum(gc_distribution) / len(gc_distribution), 2
        )
        dup_rate   = round(
            sum(1 for v in dup_tracker.values() if v > 1)
            / len(dup_tracker) * 100, 2
        ) if dup_tracker else 0
        mean_qual  = round(
            sum(
                sum(pos_quals) / len(pos_quals)
                for pos_quals in qual_per_base if pos_quals
            ) / read_length, 2
        )

        sample_results.append({
            "sample"          : sample_name,
            "reads_analyzed"  : n_reads // n_samples,
            "mean_gc_content" : mean_gc,
            "mean_quality"    : mean_qual,
            "duplication_rate": dup_rate,
            "adapter_content" : round(
                adapter_count / (n_reads // n_samples) * 100, 2
            ),
            "overrep_sequences": overrep_count,
            "qc_pass"         : (
                mean_qual >= 28 and
                mean_gc >= 40 and mean_gc <= 60 and
                dup_rate < 50
            )
        })

    elapsed = round(time.time() - start, 3)

    # MultiQC aggregation simulation
    all_pass    = sum(1 for s in sample_results if s["qc_pass"])
    mean_gc_all = round(
        sum(s["mean_gc_content"] for s in sample_results)
        / len(sample_results), 2
    )

    return {
        "step"                  : "QC_FastQC_MultiQC",
        "status"                : "✅ PASS",
        "tools"                 : ["FastQC", "MultiQC"],
        "samples_analyzed"      : n_samples,
        "reads_per_sample"      : n_reads // n_samples,
        "samples_passed_qc"     : all_pass,
        "samples_failed_qc"     : n_samples - all_pass,
        "mean_gc_content_all"   : mean_gc_all,
        "sample_details"        : sample_results,
        "multiqc_report"        : "Aggregated across all samples",
        "time_seconds"          : elapsed,
        "throughput_reads_per_sec": round(n_reads / elapsed)
    }


# ════════════════════════════════════════════════════════════
# STEP 2: TRIMMING (Trimmomatic / fastp)
# ════════════════════════════════════════════════════════════

def simulate_trimming(n_reads=100000,
                      read_length=150,
                      n_samples=6):
    """
    Simulates Trimmomatic / fastp for RNA-seq:
    - Remove adapter sequences (Illumina TruSeq)
    - Sliding window quality trimming
    - Minimum length filter
    - For RNA-seq: also remove polyA tails
    Bottleneck: CPU + SSD I/O
    """
    print("    ✂️  Step 2: Trimming (Trimmomatic / fastp)...")
    start = time.time()

    tool_results = {}

    for tool in ["Trimmomatic", "fastp"]:
        reads_passed  = 0
        reads_trimmed = 0
        reads_dropped = 0
        bases_removed = 0

        for i in range(n_reads):
            seq  = generate_rna_sequence(read_length)
            qual = generate_quality_scores(read_length)

            trimmed_len = read_length

            # Adapter trimming (TruSeq adapter)
            adapter = 'AGATCGGAAGAGC'
            if adapter in seq:
                trimmed_len = seq.index(adapter)

            # Sliding window quality trimming
            # SLIDINGWINDOW:4:15 for Trimmomatic
            # fastp uses similar approach
            window = 4
            threshold = 15 if tool == "Trimmomatic" else 20
            for j in range(0, min(trimmed_len, len(qual))-window, window):
                window_qual = sum(
                    ord(c)-33 for c in qual[j:j+window]
                ) / window
                if window_qual < threshold:
                    trimmed_len = j
                    break

            # PolyA tail removal (RNA-seq specific)
            poly_a = 'AAAAAAAAAA'
            if poly_a in seq[:trimmed_len]:
                trimmed_len = seq.index(poly_a)

            bases_removed += (read_length - trimmed_len)

            # MINLEN filter
            min_len = 36
            if trimmed_len >= min_len:
                reads_passed += 1
                if trimmed_len < read_length:
                    reads_trimmed += 1
            else:
                reads_dropped += 1

        tool_results[tool] = {
            "reads_processed" : n_reads,
            "reads_passed"    : reads_passed,
            "reads_trimmed"   : reads_trimmed,
            "reads_dropped"   : reads_dropped,
            "bases_removed"   : bases_removed,
            "survival_rate"   : round(reads_passed/n_reads*100, 2),
            "bases_removed_pct": round(
                bases_removed / (n_reads * read_length) * 100, 2
            )
        }

    elapsed = round(time.time() - start, 3)

    return {
        "step"                   : "Trimming",
        "status"                 : "✅ PASS",
        "tools"                  : ["Trimmomatic", "fastp"],
        "samples_processed"      : n_samples,
        "tool_comparison"        : tool_results,
        "time_seconds"           : elapsed,
        "throughput_reads_per_sec": round(n_reads / elapsed)
    }


# ════════════════════════════════════════════════════════════
# STEP 3: ALIGNMENT (STAR / HISAT2)
# ════════════════════════════════════════════════════════════

def simulate_star_alignment(n_reads=50000,
                             read_length=150,
                             n_samples=6):
    """
    Simulates STAR alignment for RNA-seq:
    - STAR is splice-aware (critical for RNA-seq!)
    - Detects splice junctions
    - Needs ~32GB RAM for human genome index
    - Much more RAM than BWA-MEM2
    - Outputs BAM files (~5-15GB per sample)
    Key difference from WES: splice junction detection
    """
    print("    🔗 Step 3a: Alignment Simulation (STAR)...")
    start = time.time()

    resources = get_system_resources()
    ram_gb    = resources["total_ram_gb"]

    # STAR index load simulation
    # Real STAR human genome index = ~32GB RAM
    print(f"       Simulating STAR genome index load "
          f"({ram_gb}GB RAM available, needs 32GB)...")

    index_load_start = time.time()
    star_ram_sufficient = ram_gb >= 32

    try:
        alloc_size = min(
            200_000_000,
            int(resources["available_ram_gb"] * 0.3 * 1024**3 / 8)
        )
        genome_index = np.random.randint(
            0, 4, size=alloc_size, dtype=np.uint8
        )
        index_load_time = round(time.time() - index_load_start, 3)
    except MemoryError:
        genome_index    = None
        index_load_time = None

    # Simulate splice-aware alignment
    total_aligned      = 0
    total_unmapped     = 0
    spliced_reads      = 0  # Reads spanning splice junctions
    multi_mapped       = 0
    novel_junctions    = 0

    sample_results = []

    for s in range(n_samples):
        aligned   = 0
        unmapped  = 0
        spliced   = 0
        multi     = 0
        novel_jxn = 0

        for i in range(n_reads // n_samples):
            read = generate_rna_sequence(read_length)

            if genome_index is not None:
                pos        = random.randint(
                    0, len(genome_index) - read_length
                )
                ref_chunk  = genome_index[pos:pos+read_length]
                read_array = np.frombuffer(
                    read.encode(), dtype=np.uint8
                ) % 4
                score = np.sum(read_array == ref_chunk)
            else:
                score = random.randint(0, read_length)

            # RNA-seq specific: splice junction detection
            # ~30% of RNA-seq reads span splice junctions
            is_spliced = random.random() < 0.30

            if score >= 120:
                aligned += 1
                if is_spliced:
                    spliced += 1
                    # Some junctions are novel
                    if random.random() < 0.05:
                        novel_jxn += 1
            elif score >= 90:
                multi += 1
            else:
                unmapped += 1

        sample_name = (
            f"Disease_{s+1}" if s < 3 else f"Control_{s-2}"
        )
        sample_results.append({
            "sample"           : sample_name,
            "reads_aligned"    : aligned,
            "reads_unmapped"   : unmapped,
            "multi_mapped"     : multi,
            "spliced_reads"    : spliced,
            "novel_junctions"  : novel_jxn,
            "mapping_rate"     : round(
                aligned / (n_reads // n_samples) * 100, 2
            ),
            "splicing_rate"    : round(
                spliced / aligned * 100, 2
            ) if aligned > 0 else 0
        })

        total_aligned   += aligned
        total_unmapped  += unmapped
        spliced_reads   += spliced
        multi_mapped    += multi
        novel_junctions += novel_jxn

    if genome_index is not None:
        del genome_index

    elapsed = round(time.time() - start, 3)

    return {
        "step"                   : "Alignment_STAR",
        "tool"                   : "STAR 2-pass mode",
        "star_ram_required_gb"   : 32,
        "star_ram_available_gb"  : ram_gb,
        "star_ram_sufficient"    : star_ram_sufficient,
        "index_load_time_seconds": index_load_time,
        "samples_aligned"        : n_samples,
        "total_reads_aligned"    : total_aligned,
        "total_reads_unmapped"   : total_unmapped,
        "multi_mapped"           : multi_mapped,
        "spliced_reads"          : spliced_reads,
        "novel_junctions_found"  : novel_junctions,
        "splice_aware"           : True,
        "sample_details"         : sample_results,
        "time_seconds"           : elapsed,
        "throughput_reads_per_sec": round(
            (n_reads * n_samples) / elapsed
        )
    }

def simulate_hisat2_alignment(n_reads=50000,
                               read_length=150,
                               n_samples=6):
    """
    Simulates HISAT2 alignment (alternative to STAR):
    - Also splice-aware
    - Uses only ~8GB RAM (much less than STAR)
    - Slightly less sensitive than STAR
    - Better choice for low RAM machines
    """
    print("    🔗 Step 3b: Alignment Simulation (HISAT2)...")
    start = time.time()

    sample_results = []
    total_aligned  = 0

    for s in range(n_samples):
        aligned  = 0
        unmapped = 0
        spliced  = 0
        multi    = 0

        for i in range(n_reads // n_samples):
            # HISAT2 graph-based alignment
            score      = random.gauss(128, 18)
            is_spliced = random.random() < 0.28

            if score >= 120:
                aligned += 1
                if is_spliced:
                    spliced += 1
            elif score >= 85:
                multi += 1
            else:
                unmapped += 1

        sample_name = (
            f"Disease_{s+1}" if s < 3 else f"Control_{s-2}"
        )
        sample_results.append({
            "sample"       : sample_name,
            "reads_aligned": aligned,
            "unmapped"     : unmapped,
            "multi_mapped" : multi,
            "spliced_reads": spliced,
            "mapping_rate" : round(
                aligned / (n_reads // n_samples) * 100, 2
            )
        })
        total_aligned += aligned

    elapsed = round(time.time() - start, 3)

    return {
        "step"                   : "Alignment_HISAT2",
        "tool"                   : "HISAT2",
        "ram_required_gb"        : 8,
        "splice_aware"           : True,
        "samples_aligned"        : n_samples,
        "total_reads_aligned"    : total_aligned,
        "sample_details"         : sample_results,
        "time_seconds"           : elapsed,
        "throughput_reads_per_sec": round(
            (n_reads * n_samples) / elapsed
        )
    }


# ════════════════════════════════════════════════════════════
# STEP 4: POST-PROCESSING (Samtools)
# ════════════════════════════════════════════════════════════

def simulate_post_processing(n_reads=100000, n_samples=6):
    """
    Simulates Samtools post-processing for RNA-seq:
    - Sort BAM by coordinate
    - Index BAM file
    - Remove unmapped reads
    - Get alignment statistics (flagstat)
    Note: RNA-seq post-processing is simpler than WES
    No MarkDuplicates needed (RNA-seq duplicates are biological)
    No BQSR needed
    Bottleneck: SSD I/O (sorting large BAM files)
    """
    print("    🔧 Step 4: Post-processing (Samtools)...")
    start = time.time()

    results = {}

    # ── Sort BAM ──
    sort_start = time.time()
    reads_to_sort = []
    for i in range(n_reads):
        chrom = f"chr{random.randint(1,22)}"
        pos   = random.randint(1, 250_000_000)
        reads_to_sort.append((chrom, pos))
    reads_to_sort.sort()
    results["sort_bam"] = {
        "reads_sorted"   : len(reads_to_sort),
        "time_seconds"   : round(time.time() - sort_start, 3)
    }

    # ── Index BAM ──
    index_start = time.time()
    bam_index   = {}
    for chrom, pos in reads_to_sort:
        if chrom not in bam_index:
            bam_index[chrom] = []
        bam_index[chrom].append(pos)
    results["index_bam"] = {
        "chromosomes_indexed": len(bam_index),
        "time_seconds"       : round(time.time() - index_start, 3)
    }

    # ── Remove unmapped ──
    filter_start   = time.time()
    unmapped_removed = sum(
        1 for _ in range(n_reads) if random.random() < 0.05
    )
    results["filter_unmapped"] = {
        "unmapped_removed": unmapped_removed,
        "time_seconds"    : round(time.time() - filter_start, 3)
    }

    # ── Flagstat ──
    flag_start = time.time()
    flagstat   = {
        "total_reads"    : n_reads,
        "mapped"         : n_reads - unmapped_removed,
        "properly_paired": round((n_reads - unmapped_removed) * 0.95),
        "singletons"     : round((n_reads - unmapped_removed) * 0.02),
        "mapping_rate"   : round(
            (n_reads - unmapped_removed) / n_reads * 100, 2
        )
    }
    results["flagstat"] = {
        "statistics" : flagstat,
        "time_seconds": round(time.time() - flag_start, 3)
    }

    total_time = round(sum(
        v["time_seconds"] for v in results.values()
    ), 3)
    elapsed = round(time.time() - start, 3)

    return {
        "step"                   : "Post_Processing_Samtools",
        "status"                 : "✅ PASS",
        "tool"                   : "Samtools",
        "note"                   : "No MarkDuplicates/BQSR for RNA-seq",
        "samples_processed"      : n_samples,
        "substeps"               : results,
        "flagstat"               : flagstat,
        "time_seconds"           : elapsed,
        "throughput_reads_per_sec": round(n_reads / elapsed)
    }


# ════════════════════════════════════════════════════════════
# STEP 5: QUANTIFICATION (featureCounts / Salmon)
# ════════════════════════════════════════════════════════════

def simulate_featurecounts(n_reads=100000, n_samples=6):
    """
    Simulates featureCounts quantification:
    - Count reads overlapping gene features
    - Uses GTF annotation file
    - Output: gene count matrix (genes × samples)
    - Input: BAM files from STAR/HISAT2
    Bottleneck: SSD I/O + CPU
    """
    print("    🔢 Step 5a: Quantification (featureCounts)...")
    start = time.time()

    # Simulate human genome gene list (~60,000 genes in hg38)
    n_genes      = 60000
    gene_ids     = [f"ENSG{str(i).zfill(11)}" for i in range(n_genes)]
    gene_names   = [f"GENE_{i}" for i in range(n_genes)]
    count_matrix = {}

    for s in range(n_samples):
        sample_name = (
            f"Disease_{s+1}" if s < 3 else f"Control_{s-2}"
        )
        # Generate realistic count distribution
        # Most genes have low counts, few have high counts
        counts = np.random.negative_binomial(
            1, 0.1, size=n_genes
        ).tolist()
        count_matrix[sample_name] = counts

    # Summary statistics
    all_counts     = [
        c for counts in count_matrix.values() for c in counts
    ]
    expressed_genes = sum(1 for c in all_counts if c > 0)
    mean_count      = round(sum(all_counts) / len(all_counts), 2)

    elapsed = round(time.time() - start, 3)

    return {
        "step"                   : "Quantification_featureCounts",
        "status"                 : "✅ PASS",
        "tool"                   : "featureCounts",
        "input"                  : "BAM files from STAR/HISAT2",
        "genes_in_annotation"    : n_genes,
        "expressed_genes"        : expressed_genes // n_samples,
        "mean_count_per_gene"    : mean_count,
        "samples_quantified"     : n_samples,
        "output"                 : f"Count matrix: {n_genes} genes × {n_samples} samples",
        "time_seconds"           : elapsed,
        "throughput_reads_per_sec": round(n_reads / elapsed),
        "count_matrix_shape"     : f"{n_genes} genes × {n_samples} samples"
    }

def simulate_salmon(n_reads=100000, n_samples=6):
    """
    Simulates Salmon quasi-mapping quantification:
    - Alignment-free (faster than featureCounts)
    - Works directly on FASTQ (no BAM needed)
    - Outputs TPM + estimated counts
    - Better for transcript-level quantification
    Bottleneck: CPU + RAM
    """
    print("    🔢 Step 5b: Quantification (Salmon)...")
    start = time.time()

    n_transcripts = 200000  # Human transcriptome
    sample_results = []

    for s in range(n_samples):
        sample_name = (
            f"Disease_{s+1}" if s < 3 else f"Control_{s-2}"
        )
        # Salmon outputs TPM values
        tpm_values = np.random.exponential(
            scale=10, size=n_transcripts
        )
        expressed = sum(1 for t in tpm_values if t > 0.1)

        sample_results.append({
            "sample"             : sample_name,
            "transcripts_quant"  : n_transcripts,
            "expressed_transcripts": expressed,
            "mean_tpm"           : round(float(np.mean(tpm_values)), 3),
            "mapping_rate"       : round(random.uniform(75, 92), 2)
        })

    elapsed = round(time.time() - start, 3)

    return {
        "step"                   : "Quantification_Salmon",
        "status"                 : "✅ PASS",
        "tool"                   : "Salmon",
        "method"                 : "Quasi-mapping (alignment-free)",
        "input"                  : "FASTQ files directly",
        "transcripts_quantified" : n_transcripts,
        "samples_quantified"     : n_samples,
        "output"                 : "TPM + estimated counts per transcript",
        "sample_details"         : sample_results,
        "time_seconds"           : elapsed,
        "throughput_reads_per_sec": round(n_reads / elapsed)
    }


# ════════════════════════════════════════════════════════════
# STEP 6: DIFFERENTIAL EXPRESSION (DESeq2 / edgeR)
# ════════════════════════════════════════════════════════════

def simulate_deseq2(n_genes=60000,
                    n_disease=3,
                    n_control=3):
    """
    Simulates DESeq2 differential expression analysis:
    - Normalize counts (size factors)
    - Estimate dispersions
    - Negative binomial model fitting
    - Wald test for DE
    - Multiple testing correction (BH/FDR)
    Output: DE genes with log2FC + p-value + padj
    Design: 3 disease vs 3 control samples
    Bottleneck: RAM (loading count matrix) + CPU (R)
    """
    print("    📈 Step 6a: Differential Expression (DESeq2)...")
    start = time.time()

    n_samples = n_disease + n_control
    de_genes  = []
    all_genes = []

    # Simulate count matrix
    for g in range(n_genes):
        gene_id   = f"ENSG{str(g).zfill(11)}"
        gene_name = f"GENE_{g}"

        # Base expression level
        base_expr = random.expovariate(0.1)

        # Generate counts for each sample
        disease_counts = [
            max(0, int(
                base_expr * random.gauss(1.5, 0.3)
                * random.uniform(0.5, 2.0)
            ))
            for _ in range(n_disease)
        ]
        control_counts = [
            max(0, int(
                base_expr * random.gauss(1.0, 0.3)
            ))
            for _ in range(n_control)
        ]

        # Mean counts
        mean_disease = (
            sum(disease_counts) / n_disease
            if n_disease > 0 else 0
        )
        mean_control = (
            sum(control_counts) / n_control
            if n_control > 0 else 0
        )

        # Log2 fold change
        if mean_control > 0 and mean_disease > 0:
            log2fc = np.log2(mean_disease / mean_control)
        else:
            log2fc = 0

        # Simulate Wald test p-value
        if abs(log2fc) > 1 and mean_control > 5:
            pvalue = random.uniform(0.0001, 0.05)
        else:
            pvalue = random.uniform(0.05, 1.0)

        # BH correction simulation
        padj = min(1.0, pvalue * n_genes / (g + 1))

        gene_result = {
            "gene_id"       : gene_id,
            "gene_name"     : gene_name,
            "mean_disease"  : round(mean_disease, 2),
            "mean_control"  : round(mean_control, 2),
            "log2FC"        : round(log2fc, 4),
            "pvalue"        : round(pvalue, 6),
            "padj"          : round(padj, 6),
            "significant"   : padj < 0.05 and abs(log2fc) >= 1
        }
        all_genes.append(gene_result)

        if gene_result["significant"]:
            de_genes.append(gene_result)

    # Sort by adjusted p-value
    de_genes.sort(key=lambda x: x["padj"])

    elapsed = round(time.time() - start, 3)

    up_regulated   = sum(
        1 for g in de_genes if g["log2FC"] > 0
    )
    down_regulated = sum(
        1 for g in de_genes if g["log2FC"] < 0
    )

    return {
        "step"                  : "Differential_Expression_DESeq2",
        "status"                : "✅ PASS",
        "tool"                  : "DESeq2",
        "design"                : f"{n_disease} disease vs {n_control} control",
        "total_genes_tested"    : n_genes,
        "significant_de_genes"  : len(de_genes),
        "up_regulated"          : up_regulated,
        "down_regulated"        : down_regulated,
        "top_10_de_genes"       : de_genes[:10],
        "thresholds"            : {
            "padj"   : "< 0.05",
            "log2FC" : ">= 1 or <= -1"
        },
        "time_seconds"          : elapsed,
        "throughput_reads_per_sec": round(n_genes / elapsed)
    }

def simulate_edger(n_genes=60000, n_disease=3, n_control=3):
    """
    Simulates edgeR differential expression:
    - Similar to DESeq2 but uses TMM normalization
    - Better for small sample sizes
    - Uses negative binomial + GLM
    - Quasi-likelihood F-test
    """
    print("    📈 Step 6b: Differential Expression (edgeR)...")
    start = time.time()

    de_genes = []
    for g in range(n_genes):
        base_expr = random.expovariate(0.1)
        mean_d    = base_expr * random.gauss(1.4, 0.3)
        mean_c    = base_expr * random.gauss(1.0, 0.3)

        if mean_c > 0 and mean_d > 0:
            log2fc = np.log2(abs(mean_d) / abs(mean_c))
        else:
            log2fc = 0

        pvalue = (
            random.uniform(0.0001, 0.05)
            if abs(log2fc) > 1
            else random.uniform(0.05, 1.0)
        )
        fdr = min(1.0, pvalue * n_genes / (g + 1))

        if fdr < 0.05 and abs(log2fc) >= 1:
            de_genes.append({
                "gene_id" : f"ENSG{str(g).zfill(11)}",
                "log2FC"  : round(log2fc, 4),
                "pvalue"  : round(pvalue, 6),
                "FDR"     : round(fdr, 6)
            })

    elapsed = round(time.time() - start, 3)

    return {
        "step"                  : "Differential_Expression_edgeR",
        "status"                : "✅ PASS",
        "tool"                  : "edgeR",
        "normalization"         : "TMM",
        "test"                  : "Quasi-likelihood F-test",
        "design"                : f"{n_disease} disease vs {n_control} control",
        "significant_de_genes"  : len(de_genes),
        "up_regulated"          : sum(
            1 for g in de_genes if g["log2FC"] > 0
        ),
        "down_regulated"        : sum(
            1 for g in de_genes if g["log2FC"] < 0
        ),
        "time_seconds"          : elapsed,
        "throughput_reads_per_sec": round(n_genes / elapsed)
    }


# ════════════════════════════════════════════════════════════
# STEP 7: PATHWAY ANALYSIS (clusterProfiler)
# ════════════════════════════════════════════════════════════

def simulate_pathway_analysis(n_de_genes=500):
    """
    Simulates clusterProfiler pathway analysis:
    - GO enrichment analysis (BP, MF, CC)
    - KEGG pathway enrichment
    - Fisher's exact test for enrichment
    - BH correction for multiple testing
    Output: Enriched pathways with p-values
    Bottleneck: CPU (statistical tests × pathways)
    """
    print("    🧪 Step 7: Pathway Analysis (clusterProfiler)...")
    start = time.time()

    # GO terms simulation
    go_categories = {
        "BP": "Biological Process",
        "MF": "Molecular Function",
        "CC": "Cellular Component"
    }

    # Simulate GO database
    n_go_terms  = 10000
    go_results  = {}

    for category, desc in go_categories.items():
        enriched_terms = []
        for i in range(n_go_terms // 3):
            go_id       = f"GO:{str(i).zfill(7)}"
            n_genes_term = random.randint(5, 500)
            overlap     = random.randint(
                0, min(n_de_genes, n_genes_term)
            )

            # Fisher's exact test simulation
            if overlap > 3:
                pvalue = random.uniform(0.0001, 0.1)
            else:
                pvalue = random.uniform(0.1, 1.0)

            padj = min(1.0, pvalue * n_go_terms / (i + 1))

            if padj < 0.05:
                enriched_terms.append({
                    "go_id"       : go_id,
                    "description" : f"GO_term_{category}_{i}",
                    "category"    : category,
                    "gene_count"  : overlap,
                    "bg_count"    : n_genes_term,
                    "pvalue"      : round(pvalue, 6),
                    "padj"        : round(padj, 6),
                    "rich_factor" : round(
                        overlap / n_genes_term, 4
                    )
                })

        go_results[category] = enriched_terms

    # KEGG pathway simulation
    n_kegg_pathways = 350  # Human KEGG pathways
    kegg_enriched   = []

    for i in range(n_kegg_pathways):
        pathway_id   = f"hsa{str(4000+i).zfill(5)}"
        n_genes_path = random.randint(10, 200)
        overlap      = random.randint(0, min(n_de_genes, n_genes_path))
        pvalue       = (
            random.uniform(0.0001, 0.05)
            if overlap > 5
            else random.uniform(0.05, 1.0)
        )
        padj = min(1.0, pvalue * n_kegg_pathways / (i + 1))

        if padj < 0.05:
            kegg_enriched.append({
                "pathway_id" : pathway_id,
                "description": f"KEGG_pathway_{i}",
                "gene_count" : overlap,
                "bg_count"   : n_genes_path,
                "pvalue"     : round(pvalue, 6),
                "padj"       : round(padj, 6)
            })

    elapsed = round(time.time() - start, 3)

    total_enriched = (
        sum(len(v) for v in go_results.values()) +
        len(kegg_enriched)
    )

    return {
        "step"                   : "Pathway_Analysis_clusterProfiler",
        "status"                 : "✅ PASS",
        "tool"                   : "clusterProfiler",
        "databases"              : ["GO (BP, MF, CC)", "KEGG"],
        "de_genes_input"         : n_de_genes,
        "go_terms_tested"        : n_go_terms,
        "kegg_pathways_tested"   : n_kegg_pathways,
        "enriched_go_bp"         : len(go_results.get("BP", [])),
        "enriched_go_mf"         : len(go_results.get("MF", [])),
        "enriched_go_cc"         : len(go_results.get("CC", [])),
        "enriched_kegg_pathways" : len(kegg_enriched),
        "total_enriched_terms"   : total_enriched,
        "top_kegg_pathways"      : kegg_enriched[:5],
        "time_seconds"           : elapsed,
        "throughput_reads_per_sec": round(n_go_terms / elapsed)
    }


# ════════════════════════════════════════════════════════════
# MAIN RNA-SEQ BENCHMARK
# ════════════════════════════════════════════════════════════

def run_rnaseq_benchmark(mode="3vs3"):
    """
    Run complete bulk RNA-seq pipeline benchmark
    mode: '3vs3' = 3 disease vs 3 control samples
    """
    n_samples = 6  # 3 disease + 3 control
    n_disease = 3
    n_control = 3

    print("\n🔬 BioMark — Bulk RNA-seq Pipeline Benchmark")
    print(f"   Design : {n_disease} Disease vs {n_control} Control")
    print(f"   Total  : {n_samples} samples")
    print("   Tools  : STAR/HISAT2 + featureCounts/Salmon + DESeq2/edgeR")
    print("   Note   : Splice-aware alignment (RNA-seq specific)")
    print("=" * 55)

    # ── Hardware Check ──
    resources = get_system_resources()
    warnings  = evaluate_hardware_rnaseq(resources)

    print("\n  📋 Hardware Assessment:")
    for w in warnings:
        print(f"     {w['level']} [{w['component']}] {w['message']}")
        print(f"     💡 {w['recommendation']}")
    print()

    results = {}

    # Run all pipeline steps
    results["step1_qc"]           = simulate_fastqc_multiqc(
        n_reads=60000, n_samples=n_samples
    )
    results["step2_trimming"]     = simulate_trimming(
        n_reads=60000, n_samples=n_samples
    )
    results["step3a_star"]        = simulate_star_alignment(
        n_reads=30000, n_samples=n_samples
    )
    results["step3b_hisat2"]      = simulate_hisat2_alignment(
        n_reads=30000, n_samples=n_samples
    )
    results["step4_postprocess"]  = simulate_post_processing(
        n_reads=60000, n_samples=n_samples
    )
    results["step5a_featurecounts"] = simulate_featurecounts(
        n_reads=60000, n_samples=n_samples
    )
    results["step5b_salmon"]      = simulate_salmon(
        n_reads=60000, n_samples=n_samples
    )
    results["step6a_deseq2"]      = simulate_deseq2(
        n_genes=20000,
        n_disease=n_disease,
        n_control=n_control
    )
    results["step6b_edger"]       = simulate_edger(
        n_genes=20000,
        n_disease=n_disease,
        n_control=n_control
    )
    results["step7_pathway"]      = simulate_pathway_analysis(
        n_de_genes=results["step6a_deseq2"]["significant_de_genes"]
    )

    # ── Calculate Score ──
    ram_gb = resources["total_ram_gb"]

    # RAM score — STAR needs 32GB
    ram_score   = min(100, (ram_gb / 32) * 100)

    qc_score    = min(100,
        results["step1_qc"]["throughput_reads_per_sec"] / 300
    )
    trim_score  = min(100,
        results["step2_trimming"]["throughput_reads_per_sec"] / 100000
    )
    align_score = min(100,
        results["step3a_star"]["throughput_reads_per_sec"] / 3000
    )
    quant_score = min(100,
        results["step5a_featurecounts"]["throughput_reads_per_sec"] / 30000
    )
    de_score    = min(100,
        results["step6a_deseq2"]["throughput_reads_per_sec"] / 10000
    )
    path_score  = min(100,
        results["step7_pathway"]["throughput_reads_per_sec"] / 5000
    )

    rna_score, capability, step_scores = _recalculate_honest_score(
        results, resources
    )

    total_time = round(sum(
        v.get("time_seconds", 0) for v in results.values()
    ), 2)

    # ── Final Summary ──
    print(f"\n  📊 RNA-seq Pipeline Results Summary:")
    print(f"  {'Step':<40} {'Time':>8}  {'Throughput':>15}")
    print("  " + "-" * 65)

    step_names = {
        "step1_qc"             : "Step 1: QC (FastQC + MultiQC)",
        "step2_trimming"       : "Step 2: Trimming (Trimmomatic/fastp)",
        "step3a_star"          : "Step 3a: Alignment (STAR)",
        "step3b_hisat2"        : "Step 3b: Alignment (HISAT2)",
        "step4_postprocess"    : "Step 4: Post-processing (Samtools)",
        "step5a_featurecounts" : "Step 5a: Quantification (featureCounts)",
        "step5b_salmon"        : "Step 5b: Quantification (Salmon)",
        "step6a_deseq2"        : "Step 6a: Diff. Expression (DESeq2)",
        "step6b_edger"         : "Step 6b: Diff. Expression (edgeR)",
        "step7_pathway"        : "Step 7: Pathway Analysis (clusterProfiler)",
    }

    for key, name in step_names.items():
        if key in results:
            t  = results[key].get("time_seconds", 0)
            tp = results[key].get("throughput_reads_per_sec", 0)
            print(f"  {name:<40} {t:>7.3f}s  {tp:>12,} /s")

    print("  " + "-" * 65)

    # DE Results Summary
    deseq2_de = results["step6a_deseq2"]["significant_de_genes"]
    edger_de  = results["step6b_edger"]["significant_de_genes"]
    kegg      = results["step7_pathway"]["enriched_kegg_pathways"]
    go_bp     = results["step7_pathway"]["enriched_go_bp"]

    print(f"\n  🧬 Biological Results Summary:")
    print(f"     DESeq2 significant DE genes : {deseq2_de:,}")
    print(f"     edgeR significant DE genes  : {edger_de:,}")
    print(f"     Enriched KEGG pathways      : {kegg}")
    print(f"     Enriched GO-BP terms        : {go_bp}")

    print(f"\n  ✅ RNA-seq Benchmark Complete!")
    print(f"  ⏱  Total time      : {total_time}s")
    print(f"  🧠 RAM             : {ram_gb}GB")
    print(f"  ⚡ CPU cores       : {resources['cpu_cores']}")
    print(f"  🏆 RNA-seq Score   : {rna_score}/100")
    print(f"  💡 Capability      : {capability}")

    return {
        "module"            : "RNA_seq",
        "mode"              : mode,
        "design"            : f"{n_disease}vs{n_control}",
        "score"             : rna_score,
        "capability"        : capability,
        "total_time_seconds": total_time,
        "hardware"          : resources,
        "hardware_warnings" : warnings,
        "pipeline_steps"    : results
    }


# ════════════════════════════════════════════════════════════
# HONEST SCORING PATCH
# ════════════════════════════════════════════════════════════

def _recalculate_honest_score(results, resources):
    """
    Recalculate RNA-seq score honestly:
    - If STAR fails (RAM < 32GB) → alignment step scores 0
    - If SSD insufficient → partial score
    - Score capped based on critical failures
    """
    ram_gb   = resources["total_ram_gb"]
    ssd_free = resources["free_ssd_gb"]

    star_failed = ram_gb < 32
    ssd_failed  = ssd_free < 60

    # Per-step scores
    step_scores = {
        "step1_qc"             : min(100, results["step1_qc"].get("throughput_reads_per_sec",0)/300),
        "step2_trimming"       : min(100, results["step2_trimming"].get("throughput_reads_per_sec",0)/100000),
        "step3a_star"          : 0 if star_failed else min(100, results["step3a_star"].get("throughput_reads_per_sec",0)/3000),
        "step3b_hisat2"        : min(100, results["step3b_hisat2"].get("throughput_reads_per_sec",0)/3000) if ram_gb >= 8 else 0,
        "step4_postprocess"    : min(100, results["step4_postprocess"].get("throughput_reads_per_sec",0)/50000),
        "step5a_featurecounts" : min(100, results["step5a_featurecounts"].get("throughput_reads_per_sec",0)/30000),
        "step5b_salmon"        : min(100, results["step5b_salmon"].get("throughput_reads_per_sec",0)/30000),
        "step6a_deseq2"        : min(100, results["step6a_deseq2"].get("throughput_reads_per_sec",0)/10000),
        "step6b_edger"         : min(100, results["step6b_edger"].get("throughput_reads_per_sec",0)/10000),
        "step7_pathway"        : min(100, results["step7_pathway"].get("throughput_reads_per_sec",0)/5000),
    }

    weights = {
        "step1_qc"             : 0.05,
        "step2_trimming"       : 0.05,
        "step3a_star"          : 0.20,
        "step3b_hisat2"        : 0.10,
        "step4_postprocess"    : 0.10,
        "step5a_featurecounts" : 0.10,
        "step5b_salmon"        : 0.05,
        "step6a_deseq2"        : 0.20,
        "step6b_edger"         : 0.05,
        "step7_pathway"        : 0.10,
    }

    weighted = sum(step_scores[k]*w for k,w in weights.items())

    # Cap based on failures
    if star_failed and ssd_failed:
        score      = round(min(weighted, 20), 1)
        capability = "🔴 Cannot run full RNA-seq (RAM + SSD insufficient)"
    elif star_failed:
        score      = round(min(weighted, 45), 1)
        capability = "🟡 Use HISAT2 only — STAR needs 32GB RAM"
    elif ssd_failed:
        score      = round(min(weighted, 40), 1)
        capability = "🔴 Insufficient SSD for STAR index + samples"
    elif ram_gb < 32:
        score      = round(min(weighted, 60), 1)
        capability = "🟡 HISAT2 pipeline only — no STAR"
    else:
        score      = round(weighted, 1)
        capability = "🟢 Full RNA-seq pipeline capable"

    return score, capability, step_scores
