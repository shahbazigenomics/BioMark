# BioMark - Epigenomics Benchmark Module
# Simulates ChIP-seq and ATAC-seq pipelines
# Author: Amir Shahbazi
# GitHub: shahbazigenomics

import time
import random
import os
import gzip
import multiprocessing
import psutil
import numpy as np
import warnings
warnings.filterwarnings("ignore")

# ════════════════════════════════════════════════════════════
# HELPERS
# ════════════════════════════════════════════════════════════

def get_system_resources():
    mem  = psutil.virtual_memory()
    disk = psutil.disk_usage("/")
    return {
        "total_ram_gb"    : round(mem.total / (1024**3), 1),
        "available_ram_gb": round(mem.available / (1024**3), 1),
        "total_ssd_gb"    : round(disk.total / (1024**3), 1),
        "free_ssd_gb"     : round(disk.free / (1024**3), 1),
        "cpu_cores"       : multiprocessing.cpu_count(),
    }

def evaluate_hardware_epigenomics(resources):
    warnings_list = []
    ram   = resources["total_ram_gb"]
    ssd   = resources["free_ssd_gb"]
    cores = resources["cpu_cores"]

    # RAM
    if ram < 8:
        warnings_list.append({
            "level"         : "🔴 CRITICAL", "component": "RAM",
            "message"       : f"{ram}GB — Cannot run ChIP-seq or ATAC-seq",
            "recommendation": "Minimum 16GB RAM required",
            "can_run"       : False
        })
    elif ram < 16:
        warnings_list.append({
            "level"         : "🔴 ERROR", "component": "RAM",
            "message"       : f"{ram}GB — Bowtie2 alignment will struggle",
            "recommendation": "16GB+ recommended for ChIP-seq/ATAC-seq",
            "can_run"       : False
        })
    elif ram < 32:
        warnings_list.append({
            "level"         : "🟡 WARNING", "component": "RAM",
            "message"       : f"{ram}GB — Single sample ok, multi-sample tight",
            "recommendation": "32GB for comfortable multi-sample analysis",
            "can_run"       : True
        })
    else:
        warnings_list.append({
            "level"         : "🟢 GOOD", "component": "RAM",
            "message"       : f"{ram}GB — Suitable for ChIP-seq + ATAC-seq",
            "recommendation": "Sufficient for standard epigenomics pipelines",
            "can_run"       : True
        })

    # SSD
    if ssd < 50:
        warnings_list.append({
            "level"         : "🔴 CRITICAL", "component": "SSD",
            "message"       : f"{ssd}GB free — Cannot store epigenomics data",
            "recommendation": "Need 50GB+ free per sample",
            "can_run"       : False
        })
    else:
        warnings_list.append({
            "level"         : "🟢 GOOD", "component": "SSD",
            "message"       : f"{ssd}GB free — Sufficient for epigenomics",
            "recommendation": "200GB+ for large multi-sample datasets",
            "can_run"       : True
        })

    # CPU
    ram_ok = ram >= 16
    if not ram_ok:
        warnings_list.append({
            "level"         : "⚪ N/A", "component": "CPU",
            "message"       : f"{cores} cores — Irrelevant while RAM insufficient",
            "recommendation": "Fix RAM first",
            "can_run"       : False
        })
    else:
        warnings_list.append({
            "level"         : "🟢 GOOD", "component": "CPU",
            "message"       : f"{cores} cores — Good for epigenomics",
            "recommendation": "8+ cores for faster alignment",
            "can_run"       : True
        })

    return warnings_list


# ════════════════════════════════════════════════════════════
# ChIP-seq PIPELINE
# ════════════════════════════════════════════════════════════

def simulate_chipseq_qc(n_reads=100000):
    print("    📊 ChIP Step 1: QC (FastQC)...")
    start = time.time()
    passed = sum(
        1 for _ in range(n_reads)
        if random.uniform(0,40) >= 20
    )
    elapsed = max(round(time.time()-start, 3), 0.001)
    return {
        "step"    : "ChIPseq_QC",
        "status"  : "✅ PASS",
        "reads_passed": passed,
        "time_seconds": elapsed,
        "throughput_reads_per_sec": round(n_reads/elapsed)
    }

def simulate_chipseq_alignment(resources, n_reads=100000):
    print("    🔗 ChIP Step 2: Alignment (Bowtie2)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]
    ram_ok = ram_gb >= 16  # hg38 index needs 16GB+

    if not ram_ok:
        return {
            "step"    : "ChIPseq_Alignment",
            "status"  : "🔴 FAIL",
            "note"    : f"Bowtie2 needs 16GB+ RAM for hg38. You have {ram_gb}GB.",
            "time_seconds": round(time.time()-start, 3),
            "throughput_reads_per_sec": 0
        }

    aligned = sum(1 for _ in range(n_reads) if random.random() > 0.1)
    elapsed = max(round(time.time()-start, 3), 0.001)
    return {
        "step"        : "ChIPseq_Alignment",
        "status"      : "✅ PASS",
        "reads_aligned": aligned,
        "mapping_rate": round(aligned/n_reads*100, 2),
        "time_seconds": elapsed,
        "throughput_reads_per_sec": round(n_reads/elapsed)
    }

def simulate_peak_calling(n_reads=100000):
    print("    🏔️  ChIP Step 3: Peak Calling (MACS3)...")
    start = time.time()
    chromosomes = [f"chr{i}" for i in range(1,23)] + ["chrX","chrY"]
    n_peaks     = random.randint(5000, 50000)
    peaks       = []
    for i in range(n_peaks):
        chrom  = random.choice(chromosomes)
        start_ = random.randint(1, 250_000_000)
        width  = random.randint(100, 1000)
        score  = random.uniform(10, 300)
        peaks.append({
            "chrom": chrom,
            "start": start_,
            "end"  : start_ + width,
            "score": round(score, 2)
        })
    elapsed = max(round(time.time()-start, 3), 0.001)
    return {
        "step"        : "Peak_Calling_MACS3",
        "status"      : "✅ PASS",
        "peaks_called": n_peaks,
        "narrow_peaks": int(n_peaks * 0.7),
        "broad_peaks" : int(n_peaks * 0.3),
        "mean_peak_score": round(sum(p["score"] for p in peaks)/n_peaks, 2),
        "time_seconds": elapsed,
        "throughput_reads_per_sec": round(n_peaks/elapsed)
    }

def simulate_motif_analysis(n_peaks=10000):
    print("    🔤 ChIP Step 4: Motif Analysis (HOMER)...")
    start = time.time()
    motifs = []
    for i in range(20):
        motifs.append({
            "motif_id"   : f"MOTIF_{i:03d}",
            "tf_name"    : random.choice([
                "CTCF","SP1","AP1","NF-kB","p53",
                "GATA1","REST","YY1","E2F1","RUNX1"
            ]),
            "p_value"    : round(random.uniform(1e-50, 0.05), 8),
            "pct_targets": round(random.uniform(0.1, 0.8), 3),
            "enrichment" : round(random.uniform(1.5, 20), 2)
        })
    elapsed = max(round(time.time()-start, 3), 0.001)
    return {
        "step"           : "Motif_Analysis_HOMER",
        "status"         : "✅ PASS",
        "peaks_analyzed" : n_peaks,
        "motifs_found"   : len(motifs),
        "top_motifs"     : motifs[:5],
        "time_seconds"   : elapsed,
        "throughput_reads_per_sec": round(n_peaks/elapsed)
    }

def simulate_differential_binding(n_peaks=10000):
    print("    📈 ChIP Step 5: Differential Binding (DiffBind)...")
    start = time.time()
    diff_peaks = []
    for i in range(n_peaks):
        fc   = random.gauss(0, 1.5)
        pval = random.uniform(0, 1)
        padj = min(1.0, pval * n_peaks / (i+1))
        if abs(fc) > 1 and padj < 0.05:
            diff_peaks.append({
                "peak_id": f"PEAK_{i}",
                "log2FC" : round(fc, 3),
                "padj"   : round(padj, 6)
            })
    elapsed = max(round(time.time()-start, 3), 0.001)
    return {
        "step"                  : "Differential_Binding",
        "status"                : "✅ PASS",
        "peaks_tested"          : n_peaks,
        "differential_peaks"    : len(diff_peaks),
        "gained_peaks"          : sum(1 for p in diff_peaks if p["log2FC"] > 0),
        "lost_peaks"            : sum(1 for p in diff_peaks if p["log2FC"] < 0),
        "time_seconds"          : elapsed,
        "throughput_reads_per_sec": round(n_peaks/elapsed)
    }


# ════════════════════════════════════════════════════════════
# ATAC-seq PIPELINE
# ════════════════════════════════════════════════════════════

def simulate_atacseq_qc(n_reads=100000):
    print("    📊 ATAC Step 1: QC + Fragment Size Distribution...")
    start = time.time()
    # ATAC-seq has characteristic nucleosomal pattern
    fragment_sizes = []
    for _ in range(min(n_reads, 10000)):
        # Sub-nucleosomal (<150bp), mono (~200bp), di (~400bp)
        r = random.random()
        if r < 0.5:
            fragment_sizes.append(random.randint(100, 150))
        elif r < 0.75:
            fragment_sizes.append(random.randint(180, 220))
        else:
            fragment_sizes.append(random.randint(350, 450))

    mean_frag = round(sum(fragment_sizes)/len(fragment_sizes))
    sub_nucl  = sum(1 for f in fragment_sizes if f < 150)
    elapsed   = max(round(time.time()-start, 3), 0.001)

    return {
        "step"                   : "ATACseq_QC",
        "status"                 : "✅ PASS",
        "reads_analyzed"         : n_reads,
        "mean_fragment_size_bp"  : mean_frag,
        "sub_nucleosomal_pct"    : round(sub_nucl/len(fragment_sizes)*100, 2),
        "nucleosomal_pattern"    : "✅ Detected",
        "time_seconds"           : elapsed,
        "throughput_reads_per_sec": round(n_reads/elapsed)
    }

def simulate_atacseq_alignment(resources, n_reads=100000):
    print("    🔗 ATAC Step 2: Alignment (Bowtie2)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]
    ram_ok = ram_gb >= 16  # hg38 index needs 16GB+

    if not ram_ok:
        return {
            "step"    : "ATACseq_Alignment",
            "status"  : "🔴 FAIL",
            "note"    : f"Bowtie2 needs 16GB+ RAM for hg38. You have {ram_gb}GB.",
            "time_seconds": round(time.time()-start, 3),
            "throughput_reads_per_sec": 0
        }

    aligned     = sum(1 for _ in range(n_reads) if random.random() > 0.1)
    mito_reads  = int(aligned * random.uniform(0.1, 0.4))  # Mitochondrial contamination
    usable      = aligned - mito_reads

    elapsed = max(round(time.time()-start, 3), 0.001)
    return {
        "step"              : "ATACseq_Alignment",
        "status"            : "✅ PASS",
        "reads_aligned"     : aligned,
        "mitochondrial_pct" : round(mito_reads/aligned*100, 2),
        "usable_reads"      : usable,
        "time_seconds"      : elapsed,
        "throughput_reads_per_sec": round(n_reads/elapsed)
    }

def simulate_atacseq_peaks(n_reads=100000):
    print("    🏔️  ATAC Step 3: Peak Calling (MACS3)...")
    start   = time.time()
    n_peaks = random.randint(50000, 200000)
    tss_enrichment = random.uniform(5, 20)  # TSS enrichment score
    elapsed = max(round(time.time()-start, 3), 0.001)
    return {
        "step"             : "ATACseq_Peaks",
        "status"           : "✅ PASS",
        "peaks_called"     : n_peaks,
        "tss_enrichment"   : round(tss_enrichment, 2),
        "frip_score"       : round(random.uniform(0.1, 0.6), 3),
        "time_seconds"     : elapsed,
        "throughput_reads_per_sec": round(n_peaks/elapsed)
    }

def simulate_chromvar(n_peaks=100000):
    print("    🔍 ATAC Step 4: TF Activity (chromVAR)...")
    start = time.time()
    tf_activities = []
    for i in range(300):
        tf_activities.append({
            "tf"      : f"TF_{i:03d}",
            "z_score" : round(random.gauss(0, 2), 3),
            "p_value" : round(random.uniform(0, 1), 6)
        })
    significant = sum(1 for t in tf_activities if abs(t["z_score"]) > 1.96)
    elapsed     = max(round(time.time()-start, 3), 0.001)
    return {
        "step"              : "chromVAR_TF_Activity",
        "status"            : "✅ PASS",
        "tfs_analyzed"      : len(tf_activities),
        "significant_tfs"   : significant,
        "time_seconds"      : elapsed,
        "throughput_reads_per_sec": round(len(tf_activities)/elapsed)
    }

def simulate_deeptools(n_peaks=100000):
    print("    📊 Step 5: Coverage + Visualization (deepTools)...")
    start = time.time()
    # Simulate bigWig generation and heatmap
    bins        = 500
    coverage    = np.random.poisson(10, size=bins).astype(np.float32)
    coverage    = coverage / coverage.max()
    correlation = round(random.uniform(0.7, 1.0), 3)
    elapsed     = max(round(time.time()-start, 3), 0.001)
    del coverage
    return {
        "step"              : "deepTools_Coverage",
        "status"            : "✅ PASS",
        "bins_computed"     : bins,
        "sample_correlation": correlation,
        "outputs"           : ["bigWig", "heatmap", "fingerprint plot"],
        "time_seconds"      : elapsed,
        "throughput_reads_per_sec": round(bins/elapsed)
    }


# ════════════════════════════════════════════════════════════
# MAIN EPIGENOMICS BENCHMARK
# ════════════════════════════════════════════════════════════

def run_epigenomics_benchmark():
    print("\n🧬 BioMark — Epigenomics Benchmark")
    print("   Modules : ChIP-seq + ATAC-seq")
    print("   Tools   : Bowtie2 + MACS3 + HOMER + DiffBind + chromVAR + deepTools")
    print("=" * 55)

    resources = get_system_resources()
    warnings  = evaluate_hardware_epigenomics(resources)
    ram_gb    = resources["total_ram_gb"]

    print("\n  📋 Hardware Assessment:")
    for w in warnings:
        print(f"     {w['level']} [{w['component']}] {w['message']}")
        print(f"     💡 {w['recommendation']}")
    print()

    results = {}

    # ── ChIP-seq Pipeline ──
    print("  🔬 Running ChIP-seq Pipeline...")
    results["chipseq_qc"]    = simulate_chipseq_qc()
    results["chipseq_align"] = simulate_chipseq_alignment(resources)

    chip_align_ok = "PASS" in str(results["chipseq_align"].get("status",""))
    fail_msg = {"status": "🔴 FAIL", "note": "Depends on alignment — skipped",
                "time_seconds": 0, "throughput_reads_per_sec": 0}

    results["chipseq_peaks"] = simulate_peak_calling() if chip_align_ok else {**fail_msg, "step": "Peak_Calling"}
    results["chipseq_motif"] = simulate_motif_analysis(
        results["chipseq_peaks"].get("peaks_called", 10000)
    ) if chip_align_ok else {**fail_msg, "step": "Motif_Analysis"}
    results["chipseq_diff"]  = simulate_differential_binding(
        results["chipseq_peaks"].get("peaks_called", 10000)
    ) if chip_align_ok else {**fail_msg, "step": "Differential_Binding"}

    # ── ATAC-seq Pipeline ──
    print("\n  🔬 Running ATAC-seq Pipeline...")
    results["atacseq_qc"]    = simulate_atacseq_qc()
    results["atacseq_align"] = simulate_atacseq_alignment(resources)

    atac_align_ok = "PASS" in str(results["atacseq_align"].get("status",""))

    results["atacseq_peaks"]    = simulate_atacseq_peaks() if atac_align_ok else {**fail_msg, "step": "ATAC_Peaks"}
    results["atacseq_chromvar"] = simulate_chromvar() if atac_align_ok else {**fail_msg, "step": "chromVAR"}
    results["deeptools"]        = simulate_deeptools() if (chip_align_ok or atac_align_ok) else {**fail_msg, "step": "deepTools"}

    # ── SCORING ──
    weights = {
        "chipseq_qc"      : 0.05,
        "chipseq_align"   : 0.15,
        "chipseq_peaks"   : 0.15,
        "chipseq_motif"   : 0.10,
        "chipseq_diff"    : 0.10,
        "atacseq_qc"      : 0.05,
        "atacseq_align"   : 0.15,
        "atacseq_peaks"   : 0.10,
        "atacseq_chromvar": 0.10,
        "deeptools"       : 0.05,
    }

    step_scores = {}
    for key in weights:
        s = results[key]
        if "FAIL" in str(s.get("status","")):
            step_scores[key] = 0
        else:
            tp = s.get("throughput_reads_per_sec", 0)
            step_scores[key] = min(100, tp/1000)

    weighted = sum(step_scores[k]*w for k,w in weights.items())
    fails    = sum(1 for s in results.values() if "FAIL" in str(s.get("status","")))

    if ram_gb < 16:
        score      = round(min(weighted, 20), 1)
        capability = "🔴 Cannot run ChIP-seq or ATAC-seq — RAM insufficient (needs 16GB+)"
    elif fails >= 2:
        score      = round(min(weighted, 25), 1)
        capability = "🔴 Critical steps fail — RAM or SSD insufficient"
    elif fails == 1:
        score      = round(min(weighted, 50), 1)
        capability = "🟡 Partial pipeline — some steps fail"
    elif ram_gb < 32:
        score      = round(min(weighted, 70), 1)
        capability = "🟡 Single sample analysis only"
    else:
        score      = round(weighted, 1)
        capability = "🟢 Full ChIP-seq + ATAC-seq pipeline capable"

    total_time = round(sum(
        v.get("time_seconds",0) for v in results.values()
    ), 2)

    # ── Print Summary ──
    print(f"\n  {'Step':<45} {'Status':<20} {'Time':>8}")
    print("  " + "-"*75)
    for key, r in results.items():
        name   = key.replace("_"," ").title()
        status = r.get("status","N/A")
        t      = r.get("time_seconds", 0)
        print(f"  {name:<45} {status:<20} {t:>7.3f}s")
    print("  " + "-"*75)

    # Biological summary
    chip_peaks = results["chipseq_peaks"].get("peaks_called", 0)
    diff_peaks = results["chipseq_diff"].get("differential_peaks", 0)
    atac_peaks = results["atacseq_peaks"].get("peaks_called", 0)
    sig_tfs    = results["atacseq_chromvar"].get("significant_tfs", 0)

    print(f"\n  🧬 Biological Results:")
    print(f"     ChIP-seq peaks called     : {chip_peaks:,}")
    print(f"     Differential binding peaks: {diff_peaks:,}")
    print(f"     ATAC-seq peaks called     : {atac_peaks:,}")
    print(f"     Significant TF activities : {sig_tfs}")
    print(f"\n  ✅ Epigenomics Benchmark Complete!")
    print(f"  ⏱  Total time  : {total_time}s")
    print(f"  🧠 RAM         : {ram_gb}GB")
    print(f"  🏆 Score       : {score}/100")
    print(f"  💡 Capability  : {capability}")

    return {
        "module"            : "Epigenomics",
        "score"             : score,
        "capability"        : capability,
        "total_time_seconds": total_time,
        "hardware"          : resources,
        "hardware_warnings" : warnings,
        "pipeline_steps"    : results
    }
