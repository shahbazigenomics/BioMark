# BioMark - scRNA-seq Pipeline Benchmark Module
# Simulates real single-cell RNA-seq computational workloads
# With distinct Small / Medium / Large dataset assessment
# Author: Amir Shahbazi
# GitHub: shahbazigenomics

import time
import random
import os
import multiprocessing
import psutil
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

# ════════════════════════════════════════════════════════════
# SAFE DIVISION HELPER
# ════════════════════════════════════════════════════════════

def safe_tp(numerator, elapsed):
    """Safe throughput calculation — never divide by zero"""
    return round(numerator / max(elapsed, 0.001))

# ════════════════════════════════════════════════════════════
# DATASET SIZE DEFINITIONS
# ════════════════════════════════════════════════════════════

DATASET_PROFILES = {
    "small": {
        "label"         : "Small Study",
        "n_cells"       : 5000,
        "n_genes"       : 20000,
        "n_hvg"         : 2000,
        "description"   : "~5,000 cells (e.g. single tissue biopsy)",
        "ram_minimum_gb": 16,
        "ram_optimal_gb": 32,
        "ssd_needed_gb" : 50,
        "cellranger_ram": 32,
    },
    "medium": {
        "label"         : "Medium Study",
        "n_cells"       : 30000,
        "n_genes"       : 25000,
        "n_hvg"         : 3000,
        "description"   : "~30,000 cells (e.g. PBMC or tissue atlas)",
        "ram_minimum_gb": 32,
        "ram_optimal_gb": 64,
        "ssd_needed_gb" : 200,
        "cellranger_ram": 32,
    },
    "large": {
        "label"         : "Large Study",
        "n_cells"       : 100000,
        "n_genes"       : 30000,
        "n_hvg"         : 4000,
        "description"   : "~100,000 cells (e.g. multi-sample atlas)",
        "ram_minimum_gb": 64,
        "ram_optimal_gb": 128,
        "ssd_needed_gb" : 500,
        "cellranger_ram": 64,
    },
}


# ════════════════════════════════════════════════════════════
# HELPERS
# ════════════════════════════════════════════════════════════

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


def evaluate_hardware_scrna(resources, profile):
    """
    Evaluate hardware for each dataset size separately.
    Each step has different RAM requirements.
    """
    warnings = []
    ram      = resources["total_ram_gb"]
    ssd      = resources["free_ssd_gb"]
    cores    = resources["cpu_cores"]
    label    = profile["label"]
    ram_min  = profile["ram_minimum_gb"]
    ram_opt  = profile["ram_optimal_gb"]
    ssd_need = profile["ssd_needed_gb"]
    cr_ram   = profile["cellranger_ram"]

    # ── RAM ──
    if ram < cr_ram:
        warnings.append({
            "level"         : "🔴 CRITICAL",
            "component"     : "RAM",
            "message"       : f"{ram}GB RAM — Cannot run {label} scRNA-seq",
            "recommendation": f"Cell Ranger needs {cr_ram}GB+ for {label}",
            "can_run"       : False
        })
    elif ram < ram_min:
        warnings.append({
            "level"         : "🔴 ERROR",
            "component"     : "RAM",
            "message"       : f"{ram}GB RAM — Insufficient for {label} ({ram_min}GB needed)",
            "recommendation": f"Normalization + UMAP will fail. Need {ram_min}GB+",
            "can_run"       : False
        })
    elif ram < ram_opt:
        warnings.append({
            "level"         : "🟡 WARNING",
            "component"     : "RAM",
            "message"       : f"{ram}GB RAM — Marginal for {label}",
            "recommendation": f"{ram_opt}GB recommended for comfortable {label} analysis",
            "can_run"       : True
        })
    else:
        warnings.append({
            "level"         : "🟢 GOOD",
            "component"     : "RAM",
            "message"       : f"{ram}GB RAM — Suitable for {label}",
            "recommendation": f"Sufficient for {label} scRNA-seq pipeline",
            "can_run"       : True
        })

    # ── SSD ──
    if ssd < ssd_need:
        warnings.append({
            "level"         : "🔴 CRITICAL",
            "component"     : "SSD",
            "message"       : f"{ssd}GB free — Cannot store {label} output",
            "recommendation": f"Need {ssd_need}GB+ free for {label}. Cell Ranger output alone ~50-100GB/sample",
            "can_run"       : False
        })
    else:
        warnings.append({
            "level"         : "🟢 GOOD",
            "component"     : "SSD",
            "message"       : f"{ssd}GB free — Sufficient for {label}",
            "recommendation": f"Sufficient for {label} scRNA-seq storage",
            "can_run"       : True
        })

    # ── CPU ──
    ram_ok = ram >= cr_ram
    if not ram_ok:
        warnings.append({
            "level"         : "⚪ N/A",
            "component"     : "CPU",
            "message"       : f"{cores} cores — Irrelevant while RAM is insufficient",
            "recommendation": f"Fix RAM first — need {cr_ram}GB+ before CPU matters",
            "can_run"       : False
        })
    elif cores < 8:
        warnings.append({
            "level"         : "🟡 WARNING",
            "component"     : "CPU",
            "message"       : f"{cores} cores — Cell Ranger will be slow",
            "recommendation": "8+ cores recommended",
            "can_run"       : True
        })
    elif cores < 16:
        warnings.append({
            "level"         : "🟢 GOOD",
            "component"     : "CPU",
            "message"       : f"{cores} cores — Good for {label}",
            "recommendation": "16+ cores for faster Cell Ranger",
            "can_run"       : True
        })
    else:
        warnings.append({
            "level"         : "🟢 GREAT",
            "component"     : "CPU",
            "message"       : f"{cores} cores — Excellent for {label}",
            "recommendation": "Ideal for parallel multi-sample processing",
            "can_run"       : True
        })

    return warnings


# ════════════════════════════════════════════════════════════
# STEP 1: CELL RANGER
# ════════════════════════════════════════════════════════════

def simulate_cellranger(resources, profile):
    """
    Cell Ranger alignment + cell calling
    RAM: 32GB (small/medium), 64GB (large)
    SSD: 50-100GB output per sample
    """
    print("    🔬 Step 1: Cell Ranger (alignment + cell calling)...")
    start    = time.time()
    ram_gb   = resources["total_ram_gb"]
    n_reads  = profile["n_cells"] * 10
    n_cells  = profile["n_cells"]
    cr_ram   = profile["cellranger_ram"]
    ram_ok   = ram_gb >= cr_ram

    if not ram_ok:
        return {
            "step"        : "Cell_Ranger",
            "status"      : "🔴 FAIL",
            "note"        : f"Cell Ranger needs {cr_ram}GB RAM for {profile['label']}. You have {ram_gb}GB.",
            "ram_required": f"{cr_ram}GB",
            "ram_available": f"{ram_gb}GB",
            "time_seconds": round(time.time()-start, 3),
            "throughput_reads_per_sec": 0
        }

    aligned      = sum(1 for _ in range(min(n_reads, 50000)) if random.random() > 0.15)
    cells_called = int(n_cells * random.uniform(0.90, 1.10))
    median_genes = random.randint(1500, 4000)
    median_umis  = random.randint(3000, 15000)

    elapsed = max(round(time.time()-start, 3), 0.001)
    tp      = safe_tp(min(n_reads, 50000), elapsed)

    return {
        "step"                  : "Cell_Ranger",
        "status"                : "✅ PASS",
        "reads_processed"       : n_reads,
        "cells_called"          : cells_called,
        "median_genes_per_cell" : median_genes,
        "median_umis_per_cell"  : median_umis,
        "mapping_rate"          : round(aligned/min(n_reads,50000)*100, 2),
        "ram_required"          : f"{cr_ram}GB",
        "ram_available"         : f"{ram_gb}GB",
        "time_seconds"          : elapsed,
        "throughput_reads_per_sec": tp
    }


# ════════════════════════════════════════════════════════════
# STEP 2: QC
# ════════════════════════════════════════════════════════════

def simulate_scrna_qc(resources, profile):
    """
    scRNA-seq QC filtering
    RAM needed: ~2-4GB (small), ~8-16GB (medium), ~32-64GB (large)
    """
    print("    📊 Step 2: QC + Cell Filtering (Seurat/Scanpy)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]
    n_cells = profile["n_cells"]

    # RAM needed per dataset size
    ram_needed = {
        "small" : 4,
        "medium": 16,
        "large" : 48
    }
    size_key = next(
        k for k,v in DATASET_PROFILES.items()
        if v["label"] == profile["label"]
    )
    ram_for_qc = ram_needed[size_key]
    ram_ok     = ram_gb >= ram_for_qc

    if not ram_ok:
        return {
            "step"        : "scRNA_QC",
            "status"      : "🔴 FAIL",
            "note"        : f"QC needs {ram_for_qc}GB RAM for {profile['label']}. You have {ram_gb}GB.",
            "time_seconds": round(time.time()-start, 3),
            "throughput_reads_per_sec": 0
        }

    passed = dropped = low_gene = high_mito = doublets = 0
    for i in range(min(n_cells, 10000)):
        n_g      = random.randint(200, 8000)
        pct_mito = random.uniform(0, 50)
        n_umi    = random.randint(500, 50000)
        is_dbl   = random.random() < 0.05

        if is_dbl:
            doublets += 1; dropped += 1
        elif not (200 <= n_g <= 6000):
            low_gene += 1; dropped += 1
        elif pct_mito >= 20:
            high_mito += 1; dropped += 1
        elif n_umi < 500:
            dropped += 1
        else:
            passed += 1

    elapsed = max(round(time.time()-start, 3), 0.001)
    tp      = safe_tp(min(n_cells, 10000), elapsed)
    surv    = round(passed/min(n_cells,10000)*100, 2)

    return {
        "step"            : "scRNA_QC",
        "status"          : "✅ PASS",
        "cells_input"     : n_cells,
        "cells_passed"    : int(n_cells * surv/100),
        "survival_rate"   : surv,
        "removed_doublets": doublets,
        "removed_low_genes": low_gene,
        "removed_high_mito": high_mito,
        "ram_needed_gb"   : ram_for_qc,
        "time_seconds"    : elapsed,
        "throughput_reads_per_sec": tp
    }


# ════════════════════════════════════════════════════════════
# STEP 3: NORMALIZATION + HVG
# ════════════════════════════════════════════════════════════

def simulate_normalization(resources, profile, cells_after_qc):
    """
    Normalization + Highly Variable Gene selection
    RAM: full cell×gene matrix in memory
    Small:  5k  × 20k genes × 4 bytes = 400MB  → needs 4GB
    Medium: 30k × 25k genes × 4 bytes = 3GB    → needs 16GB
    Large:  100k× 30k genes × 4 bytes = 12GB   → needs 64GB
    """
    print("    🔢 Step 3: Normalization + HVG Selection...")
    start    = time.time()
    ram_gb   = resources["total_ram_gb"]
    n_genes  = profile["n_genes"]
    n_hvg    = profile["n_hvg"]

    # Matrix size calculation
    matrix_gb  = (cells_after_qc * n_genes * 4) / (1024**3)
    # Need ~5× matrix size for normalization intermediates
    ram_needed = round(matrix_gb * 5, 1)
    ram_ok     = ram_gb >= ram_needed

    if not ram_ok:
        return {
            "step"          : "Normalization",
            "status"        : "🔴 FAIL",
            "note"          : f"Matrix ({cells_after_qc:,} cells × {n_genes:,} genes) needs ~{ram_needed}GB RAM. You have {ram_gb}GB.",
            "matrix_size_gb": round(matrix_gb, 2),
            "ram_needed_gb" : ram_needed,
            "time_seconds"  : round(time.time()-start, 3),
            "throughput_reads_per_sec": 0
        }

    # Simulate normalization on reduced size
    sim_cells = min(cells_after_qc, 2000)
    sim_genes = min(n_genes, 2000)
    counts    = np.random.negative_binomial(
        1, 0.9, size=(sim_cells, sim_genes)
    ).astype(np.float32)

    # LogNormalize
    row_sums = counts.sum(axis=1, keepdims=True)
    row_sums[row_sums==0] = 1
    norm     = np.log1p(counts / row_sums * 10000)

    # HVG selection — guard against empty arrays
    if norm.shape[0] > 0 and norm.shape[1] > 0:
        means       = np.nanmean(norm, axis=0)
        variances   = np.nanvar(norm, axis=0)
        disp        = variances / (means + 1e-9)
        disp        = np.nan_to_num(disp, nan=0.0)
        n_hvg_found = int(np.sum(disp > np.percentile(disp, 75)))
    else:
        n_hvg_found = 0

    del counts, norm
    elapsed = max(round(time.time()-start, 3), 0.001)
    tp      = safe_tp(cells_after_qc, elapsed)

    return {
        "step"               : "Normalization",
        "status"             : "✅ PASS",
        "cells_normalized"   : cells_after_qc,
        "genes_tested"       : n_genes,
        "hvg_selected"       : n_hvg,
        "method"             : "LogNormalize + ScaleData",
        "matrix_size_gb"     : round(matrix_gb, 2),
        "ram_needed_gb"      : ram_needed,
        "time_seconds"       : elapsed,
        "throughput_reads_per_sec": tp
    }


# ════════════════════════════════════════════════════════════
# STEP 4: PCA + UMAP
# ════════════════════════════════════════════════════════════

def simulate_dim_reduction(resources, profile, cells_after_qc):
    """
    PCA + UMAP dimensionality reduction
    RAM: PCA matrix = cells × HVG × 4 bytes
    Small:  5k  × 2k HVG = 40MB  → needs 2GB
    Medium: 30k × 3k HVG = 360MB → needs 8GB
    Large:  100k× 4k HVG = 1.6GB → needs 32GB
    UMAP on large datasets: very slow + RAM intensive
    """
    print("    📉 Step 4: PCA + UMAP (dimensionality reduction)...")
    start   = time.time()
    ram_gb  = resources["total_ram_gb"]
    n_hvg   = profile["n_hvg"]

    # RAM needed for PCA + UMAP
    pca_matrix_gb = (cells_after_qc * n_hvg * 4) / (1024**3)
    ram_needed    = round(max(2, pca_matrix_gb * 8), 1)
    ram_ok        = ram_gb >= ram_needed

    if not ram_ok:
        return {
            "step"        : "PCA_UMAP",
            "status"      : "🔴 FAIL",
            "note"        : f"PCA needs ~{ram_needed}GB RAM for {profile['label']}. You have {ram_gb}GB.",
            "ram_needed_gb": ram_needed,
            "time_seconds": round(time.time()-start, 3),
            "throughput_reads_per_sec": 0
        }

    # Simulate PCA — guard for empty cells
    sim_size = max(min(cells_after_qc, 1000), 2)
    sim_hvg  = max(min(n_hvg, 500), 2)
    data     = np.random.randn(sim_size, sim_hvg).astype(np.float32)
    U, S, Vt = np.linalg.svd(data, full_matrices=False)
    n_pcs    = min(50, S.shape[0])
    del data, U, S, Vt

    # UMAP time estimate
    umap_note = ""
    if cells_after_qc > 50000:
        umap_note = f"⚠️ UMAP on {cells_after_qc:,} cells will take 2-6 hours on real data"
    elif cells_after_qc > 20000:
        umap_note = f"UMAP on {cells_after_qc:,} cells: ~30-60 min on real data"
    else:
        umap_note = f"UMAP on {cells_after_qc:,} cells: ~5-15 min on real data"

    elapsed = max(round(time.time()-start, 3), 0.001)
    tp      = safe_tp(cells_after_qc, elapsed)

    return {
        "step"            : "PCA_UMAP",
        "status"          : "✅ PASS",
        "cells_processed" : cells_after_qc,
        "n_pcs"           : n_pcs,
        "n_hvg_input"     : n_hvg,
        "ram_needed_gb"   : ram_needed,
        "umap_note"       : umap_note,
        "time_seconds"    : elapsed,
        "throughput_reads_per_sec": tp
    }


# ════════════════════════════════════════════════════════════
# STEP 5: CLUSTERING
# ════════════════════════════════════════════════════════════

def simulate_clustering(resources, profile, cells_after_qc):
    """
    Leiden/Louvain clustering
    RAM: KNN graph = cells × neighbors × 4 bytes
    Small:  5k  × 15 neighbors = 300KB  → fine
    Medium: 30k × 15 neighbors = 1.8MB  → fine
    Large:  100k× 20 neighbors = 8MB    → fine
    Not RAM-limited but CPU-limited for large datasets
    """
    print("    🔵 Step 5: Clustering (Leiden algorithm)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]

    # Clustering is not RAM limited
    # but large datasets need more CPU time
    # Guard: if no cells passed QC, return FAIL
    if cells_after_qc == 0:
        return {
            "step"      : "Clustering",
            "status"    : "🔴 FAIL",
            "note"      : "No cells available after QC/normalization",
            "n_clusters": 0,
            "time_seconds": round(time.time()-start, 3),
            "throughput_reads_per_sec": 0
        }

    n_clusters  = random.randint(8, 25)
    assignments = [random.randint(0, n_clusters-1)
                   for _ in range(min(cells_after_qc, 5000))]
    sizes       = {}
    for c in assignments:
        sizes[c] = sizes.get(c, 0) + 1

    # Guard against empty sizes
    if not sizes:
        sizes = {0: cells_after_qc}

    elapsed = max(round(time.time()-start, 3), 0.001)
    tp      = safe_tp(min(cells_after_qc, 5000), elapsed)

    # Real time estimate
    if cells_after_qc > 50000:
        time_note = f"Real clustering of {cells_after_qc:,} cells: ~30-60 min"
    elif cells_after_qc > 20000:
        time_note = f"Real clustering of {cells_after_qc:,} cells: ~5-15 min"
    else:
        time_note = f"Real clustering of {cells_after_qc:,} cells: ~1-3 min"

    return {
        "step"            : "Clustering",
        "status"          : "✅ PASS",
        "cells_clustered" : cells_after_qc,
        "n_clusters"      : n_clusters,
        "algorithm"       : "Leiden",
        "resolution"      : 0.5,
        "largest_cluster" : max(sizes.values()),
        "smallest_cluster": min(sizes.values()),
        "time_note"       : time_note,
        "time_seconds"    : elapsed,
        "throughput_reads_per_sec": tp
    }


# ════════════════════════════════════════════════════════════
# STEP 6: MARKER GENES
# ════════════════════════════════════════════════════════════

def simulate_marker_genes(resources, profile, cells_after_qc, n_clusters):
    """
    Marker gene identification (Wilcoxon test)
    RAM: needs full expression matrix per cluster comparison
    Small:  manageable with 16GB
    Medium: needs 32GB+
    Large:  needs 64GB+ (100k cells × 30k genes)
    """
    print("    🏷️  Step 6: Marker Genes + Cell Type Annotation...")
    start    = time.time()
    ram_gb   = resources["total_ram_gb"]
    n_genes  = profile["n_genes"]

    # RAM needed: matrix + comparison intermediates
    matrix_gb  = (cells_after_qc * n_genes * 4) / (1024**3)
    ram_needed = round(matrix_gb * 3, 1)
    ram_ok     = ram_gb >= ram_needed

    if not ram_ok:
        return {
            "step"        : "Marker_Genes",
            "status"      : "🔴 FAIL",
            "note"        : f"Wilcoxon tests need ~{ram_needed}GB RAM. You have {ram_gb}GB.",
            "ram_needed_gb": ram_needed,
            "time_seconds": round(time.time()-start, 3),
            "throughput_reads_per_sec": 0
        }

    cell_types = [
        "T cells", "B cells", "NK cells", "Monocytes",
        "Dendritic cells", "Macrophages", "Neutrophils",
        "Plasma cells", "Fibroblasts", "Endothelial cells",
        "Epithelial cells", "Stem cells", "Mast cells",
        "Basophils", "Eosinophils"
    ]

    if cells_after_qc == 0 or n_clusters == 0:
        return {
            "step"      : "Marker_Genes",
            "status"    : "🔴 FAIL",
            "note"      : "No cells or clusters available for marker gene analysis",
            "time_seconds": round(time.time()-start, 3),
            "throughput_reads_per_sec": 0
        }

    markers = {}
    total_markers = 0
    for c in range(n_clusters):
        n_markers = random.randint(50, 300)
        total_markers += n_markers
        markers[f"Cluster_{c}"] = {
            "n_markers"  : n_markers,
            "cell_type"  : cell_types[c % len(cell_types)],
            "top_markers": [f"GENE_{random.randint(1,n_genes)}"
                           for _ in range(5)]
        }

    elapsed = max(round(time.time()-start, 3), 0.001)
    tp      = safe_tp(n_genes * n_clusters, elapsed)

    return {
        "step"              : "Marker_Genes",
        "status"            : "✅ PASS",
        "clusters_analyzed" : n_clusters,
        "genes_tested"      : n_genes,
        "total_markers"     : total_markers,
        "ram_needed_gb"     : ram_needed,
        "cluster_annotations": markers,
        "time_seconds"      : elapsed,
        "throughput_reads_per_sec": tp
    }


# ════════════════════════════════════════════════════════════
# STEP 7: TRAJECTORY (Monocle3)
# ════════════════════════════════════════════════════════════

def simulate_trajectory(resources, profile, cells_after_qc):
    """
    Monocle3 trajectory analysis
    RAM: similar to clustering — not heavily RAM limited
    but large datasets need significant CPU time
    """
    print("    🛤️  Step 7: Trajectory Analysis (Monocle3)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]

    if cells_after_qc == 0:
        return {
            "step"      : "Trajectory_Monocle3",
            "status"    : "🔴 FAIL",
            "note"      : "No cells available for trajectory analysis",
            "time_seconds": round(time.time()-start, 3),
            "throughput_reads_per_sec": 0
        }

    n_nodes      = random.randint(50, 200)
    pseudotimes  = [random.uniform(0, 20)
                    for _ in range(min(cells_after_qc, 5000))]
    traj_genes   = random.randint(200, 1000)

    # Real time estimate
    if cells_after_qc > 50000:
        time_note = f"Real trajectory of {cells_after_qc:,} cells: ~2-4 hours"
    else:
        time_note = f"Real trajectory of {cells_after_qc:,} cells: ~15-30 min"

    elapsed = max(round(time.time()-start, 3), 0.001)
    tp      = safe_tp(min(cells_after_qc, 5000), elapsed)

    return {
        "step"              : "Trajectory_Monocle3",
        "status"            : "✅ PASS",
        "cells_ordered"     : cells_after_qc,
        "n_trajectory_nodes": n_nodes,
        "trajectory_genes"  : traj_genes,
        "pseudotime_range"  : "0 → 20",
        "time_note"         : time_note,
        "time_seconds"      : elapsed,
        "throughput_reads_per_sec": tp
    }


# ════════════════════════════════════════════════════════════
# MAIN scRNA-seq BENCHMARK
# ════════════════════════════════════════════════════════════

def run_scrna_benchmark():
    """
    Run scRNA-seq benchmark for all three dataset sizes:
    Small (5k cells), Medium (30k cells), Large (100k cells)
    """
    print("\n🔬 BioMark — scRNA-seq Pipeline Benchmark")
    print("   Platform : 10x Genomics Chromium")
    print("   Tools    : Cell Ranger + Seurat/Scanpy + Monocle3")
    print("   Sizes    : Small (5k) | Medium (30k) | Large (100k)")
    print("=" * 60)

    resources = get_system_resources()
    ram_gb    = resources["total_ram_gb"]
    ssd_free  = resources["free_ssd_gb"]
    cores     = resources["cpu_cores"]

    all_size_results = {}

    for size_key, profile in DATASET_PROFILES.items():

        print(f"\n{'='*60}")
        print(f"  📊 {profile['label']} — {profile['description']}")
        print(f"{'='*60}")

        # Hardware assessment for this size
        warnings = evaluate_hardware_scrna(resources, profile)
        print(f"\n  📋 Hardware Assessment for {profile['label']}:")
        for w in warnings:
            print(f"     {w['level']} [{w['component']}] {w['message']}")
            print(f"     💡 {w['recommendation']}")
        print()

        results = {}

        # Step 1: Cell Ranger
        results["step1_cellranger"] = simulate_cellranger(
            resources, profile
        )

        # Get cells after QC
        cr_ok = results["step1_cellranger"]["status"] == "✅ PASS"
        n_cells_qc = int(profile["n_cells"] * 0.80) if cr_ok else 0

        # Step 2: QC
        results["step2_qc"] = simulate_scrna_qc(resources, profile)
        if results["step2_qc"]["status"] == "✅ PASS":
            n_cells_qc = results["step2_qc"].get(
                "cells_passed", n_cells_qc
            )

        # Step 3: Normalization
        results["step3_norm"] = simulate_normalization(
            resources, profile, n_cells_qc
        )

        # Step 4: PCA + UMAP
        results["step4_dimred"] = simulate_dim_reduction(
            resources, profile, n_cells_qc
        )

        # Step 5: Clustering
        results["step5_clustering"] = simulate_clustering(
            resources, profile, n_cells_qc
        )
        n_clusters = results["step5_clustering"]["n_clusters"]

        # Step 6: Marker genes
        results["step6_markers"] = simulate_marker_genes(
            resources, profile, n_cells_qc, n_clusters
        )

        # Step 7: Trajectory
        results["step7_trajectory"] = simulate_trajectory(
            resources, profile, n_cells_qc
        )

        # ── HONEST SCORING ──
        step_scores = {}
        weights = {
            "step1_cellranger" : 0.30,
            "step2_qc"         : 0.10,
            "step3_norm"       : 0.20,
            "step4_dimred"     : 0.15,
            "step5_clustering" : 0.10,
            "step6_markers"    : 0.10,
            "step7_trajectory" : 0.05,
        }

        for key in weights:
            s = results[key]
            if "FAIL" in str(s.get("status","")):
                step_scores[key] = 0
            else:
                tp = s.get("throughput_reads_per_sec", 0)
                step_scores[key] = min(100, tp/1000)

        weighted = sum(
            step_scores[k]*w for k,w in weights.items()
        )

        # Cap based on failures and RAM
        fails = sum(
            1 for s in results.values()
            if "FAIL" in str(s.get("status",""))
        )

        if fails >= 3:
            score      = round(min(weighted, 10), 1)
            capability = f"🔴 Cannot run {profile['label']} scRNA-seq"
        elif fails >= 1:
            score      = round(min(weighted, 30), 1)
            capability = f"🔴 Critical steps fail for {profile['label']}"
        elif ram_gb < profile["ram_optimal_gb"]:
            score      = round(min(weighted, 65), 1)
            capability = f"🟡 Marginal for {profile['label']} — upgrade RAM recommended"
        else:
            score      = round(weighted, 1)
            capability = f"🟢 Can run full {profile['label']} scRNA-seq pipeline"

        total_time = round(sum(
            v.get("time_seconds",0) for v in results.values()
        ), 2)

        # ── Print Step Summary ──
        step_names = {
            "step1_cellranger" : "Step 1: Cell Ranger",
            "step2_qc"         : "Step 2: QC Filtering",
            "step3_norm"       : "Step 3: Normalization + HVG",
            "step4_dimred"     : "Step 4: PCA + UMAP",
            "step5_clustering" : "Step 5: Leiden Clustering",
            "step6_markers"    : "Step 6: Marker Genes",
            "step7_trajectory" : "Step 7: Trajectory (Monocle3)",
        }

        print(f"  {'Step':<35} {'Status':<25} {'Time':>8}")
        print("  " + "-"*70)
        for key, name in step_names.items():
            r      = results[key]
            status = r.get("status","N/A")
            t      = r.get("time_seconds",0)
            note   = r.get("note","")
            note_str = f" — {note[:50]}" if note else ""
            print(f"  {name:<35} {status:<25} {t:>7.3f}s{note_str}")
        print("  " + "-"*70)
        print(f"  💡 {capability}")
        print(f"  🏆 Score: {score}/100  |  Time: {total_time}s")

        all_size_results[size_key] = {
            "profile"           : profile,
            "score"             : score,
            "capability"        : capability,
            "total_time_seconds": total_time,
            "hardware_warnings" : warnings,
            "pipeline_steps"    : results
        }

    # ── Overall Summary ──
    print(f"\n{'='*60}")
    print(f"  📊 scRNA-seq Summary — {ram_gb}GB RAM | {ssd_free}GB SSD free")
    print(f"{'='*60}")
    for size_key, data in all_size_results.items():
        p     = data["profile"]
        score = data["score"]
        cap   = data["capability"]
        print(f"  {p['label']:<20} Score: {score:>5}/100  {cap}")
    print(f"{'='*60}")

    # Overall score = average of three sizes
    overall_score = round(
        sum(d["score"] for d in all_size_results.values()) / 3, 1
    )
    overall_cap = all_size_results["small"]["capability"]

    return {
        "module"            : "scRNA_seq",
        "score"             : overall_score,
        "capability"        : overall_cap,
        "total_time_seconds": round(sum(
            d["total_time_seconds"]
            for d in all_size_results.values()
        ), 2),
        "hardware"          : resources,
        "size_results"      : all_size_results,
        "pipeline_steps"    : all_size_results["small"]["pipeline_steps"]
    }
