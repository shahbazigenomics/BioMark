# BioMark - Long Read Sequencing Benchmark
# Simulates Nanopore + PacBio pipelines
# Author: Amir Shahbazi
# GitHub: shahbazigenomics

import time
import random
import multiprocessing
import psutil
import numpy as np
import warnings
warnings.filterwarnings("ignore")

def get_system_resources():
    mem  = psutil.virtual_memory()
    disk = psutil.disk_usage("/")
    return {
        "total_ram_gb"    : round(mem.total / (1024**3), 1),
        "available_ram_gb": round(mem.available / (1024**3), 1),
        "free_ssd_gb"     : round(disk.free / (1024**3), 1),
        "cpu_cores"       : multiprocessing.cpu_count(),
    }

def safe_tp(n, elapsed):
    return round(n / max(elapsed, 0.001))

def evaluate_hardware(resources):
    warnings_list = []
    ram  = resources["total_ram_gb"]
    ssd  = resources["free_ssd_gb"]

    if ram < 16:
        warnings_list.append({"level":"🔴 CRITICAL","component":"RAM",
            "message":f"{ram}GB — Cannot run Minimap2 for long reads",
            "recommendation":"Minimum 16GB RAM for long read alignment","can_run":False})
    elif ram < 32:
        warnings_list.append({"level":"🟡 WARNING","component":"RAM",
            "message":f"{ram}GB — Long read assembly will be tight",
            "recommendation":"32GB+ for genome assembly with Flye","can_run":True})
    else:
        warnings_list.append({"level":"🟢 GOOD","component":"RAM",
            "message":f"{ram}GB — Suitable for long read analysis",
            "recommendation":"64GB+ for large genome assembly","can_run":True})

    if ssd < 100:
        warnings_list.append({"level":"🔴 CRITICAL","component":"SSD",
            "message":f"{ssd}GB free — Cannot store long read data",
            "recommendation":"Long read FASTQ files are very large (20-100GB per run)","can_run":False})
    else:
        warnings_list.append({"level":"🟢 GOOD","component":"SSD",
            "message":f"{ssd}GB free — Sufficient","recommendation":"Sufficient","can_run":True})

    ram_ok = ram >= 16
    warnings_list.append({
        "level": "⚪ N/A" if not ram_ok else "🟢 GOOD",
        "component": "CPU",
        "message": f"{resources['cpu_cores']} cores — {'Irrelevant while RAM insufficient' if not ram_ok else 'Good for long read analysis'}",
        "recommendation": "Fix RAM first" if not ram_ok else "16+ cores for faster Minimap2",
        "can_run": ram_ok
    })
    return warnings_list

def simulate_basecalling(n_reads=10000):
    print("    📡 Step 1: Basecalling (Guppy/Dorado)...")
    start = time.time()
    # Simulate Oxford Nanopore basecalling
    reads = []
    for i in range(min(n_reads, 5000)):
        length   = random.randint(500, 50000)  # Long reads!
        quality  = random.uniform(7, 15)       # Q7-Q15 typical for Nanopore
        reads.append({"length":length, "quality":round(quality,1)})

    mean_len   = round(sum(r["length"] for r in reads)/len(reads))
    mean_qual  = round(sum(r["quality"] for r in reads)/len(reads), 1)
    n50        = sorted([r["length"] for r in reads])[int(len(reads)*0.5)]
    elapsed    = max(round(time.time()-start,3), 0.001)

    return {"step":"Basecalling","status":"✅ PASS",
            "reads_basecalled":n_reads,"mean_read_length_bp":mean_len,
            "mean_quality":mean_qual,"n50_bp":n50,
            "tool":"Dorado (Oxford Nanopore)",
            "time_seconds":elapsed,"throughput_reads_per_sec":safe_tp(n_reads,elapsed)}

def simulate_qc(n_reads=10000):
    print("    📊 Step 2: QC (NanoPlot + FastQC)...")
    start = time.time()
    passed = int(n_reads * random.uniform(0.85, 0.95))
    elapsed = max(round(time.time()-start,3), 0.001)
    return {"step":"LongRead_QC","status":"✅ PASS",
            "reads_passed":passed,"survival_rate":round(passed/n_reads*100,2),
            "min_length_filter":"500bp","min_quality_filter":"Q7",
            "time_seconds":elapsed,"throughput_reads_per_sec":safe_tp(n_reads,elapsed)}

def simulate_minimap2(resources, n_reads=10000):
    print("    🔗 Step 3: Alignment (Minimap2)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]
    ram_ok = ram_gb >= 16

    if not ram_ok:
        return {"step":"Minimap2","status":"🔴 FAIL",
                "note":f"Minimap2 needs 16GB+ RAM for hg38. You have {ram_gb}GB.",
                "time_seconds":round(time.time()-start,3),"throughput_reads_per_sec":0}

    aligned = int(n_reads * random.uniform(0.85, 0.98))
    elapsed = max(round(time.time()-start,3), 0.001)
    return {"step":"Minimap2","status":"✅ PASS",
            "reads_aligned":aligned,"mapping_rate":round(aligned/n_reads*100,2),
            "preset":"map-ont (Nanopore) / map-pb (PacBio)",
            "time_seconds":elapsed,"throughput_reads_per_sec":safe_tp(n_reads,elapsed)}

def simulate_sv_calling(resources, n_reads=10000):
    print("    🧬 Step 4: Structural Variant Calling (Sniffles2)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]
    if ram_gb < 16:
        return {"step":"SV_Calling","status":"🔴 FAIL",
                "note":"Depends on alignment — skipped",
                "time_seconds":0,"throughput_reads_per_sec":0}

    sv_types = ["DEL","INS","DUP","INV","TRA"]
    n_svs    = random.randint(5000, 30000)
    svs      = [{"type":random.choice(sv_types),
                 "size":random.randint(50,1000000),
                 "qual":random.uniform(10,60)} for _ in range(n_svs)]
    elapsed = max(round(time.time()-start,3), 0.001)
    return {"step":"SV_Calling","status":"✅ PASS",
            "svs_called":n_svs,
            "deletions":sum(1 for s in svs if s["type"]=="DEL"),
            "insertions":sum(1 for s in svs if s["type"]=="INS"),
            "duplications":sum(1 for s in svs if s["type"]=="DUP"),
            "tool":"Sniffles2","time_seconds":elapsed,
            "throughput_reads_per_sec":safe_tp(n_svs,elapsed)}

def simulate_medaka(resources, n_reads=10000):
    print("    🔧 Step 5: Variant Calling + Polishing (Medaka)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]
    if ram_gb < 16:
        return {"step":"Medaka","status":"🔴 FAIL",
                "note":"Depends on alignment — skipped",
                "time_seconds":0,"throughput_reads_per_sec":0}

    variants   = random.randint(1000, 50000)
    snps       = int(variants * 0.7)
    indels     = variants - snps
    elapsed    = max(round(time.time()-start,3), 0.001)
    return {"step":"Medaka","status":"✅ PASS",
            "variants_called":variants,"snps":snps,"indels":indels,
            "consensus_accuracy":"~99.5% (Q23)",
            "time_seconds":elapsed,"throughput_reads_per_sec":safe_tp(variants,elapsed)}

def simulate_methylation(resources):
    print("    🧪 Step 6: Methylation Detection (Modbam2bed)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]
    if ram_gb < 16:
        return {"step":"Methylation","status":"🔴 FAIL",
                "note":"Depends on alignment — skipped",
                "time_seconds":0,"throughput_reads_per_sec":0}

    # Nanopore specific — can detect 5mC methylation directly
    n_sites       = random.randint(1000000, 5000000)
    methylated    = int(n_sites * random.uniform(0.3, 0.7))
    elapsed       = max(round(time.time()-start,3), 0.001)
    return {"step":"Methylation","status":"✅ PASS",
            "cpg_sites_analyzed":n_sites,"methylated_sites":methylated,
            "methylation_rate":round(methylated/n_sites*100,2),
            "modification":"5mC (CpG methylation)",
            "note":"Unique to long read sequencing — not possible with short reads!",
            "time_seconds":elapsed,"throughput_reads_per_sec":safe_tp(n_sites,elapsed)}

def run_longread_benchmark():
    print("\n🧬 BioMark — Long Read Sequencing Benchmark")
    print("   Platforms: Oxford Nanopore + PacBio HiFi")
    print("   Tools    : Dorado + Minimap2 + Sniffles2 + Medaka")
    print("=" * 55)

    resources = get_system_resources()
    warnings  = evaluate_hardware(resources)
    ram_gb    = resources["total_ram_gb"]

    print("\n  📋 Hardware Assessment:")
    for w in warnings:
        print(f"     {w['level']} [{w['component']}] {w['message']}")
        print(f"     💡 {w['recommendation']}")
    print()

    results = {}
    results["step1_basecalling"] = simulate_basecalling()
    results["step2_qc"]          = simulate_qc()
    results["step3_minimap2"]    = simulate_minimap2(resources)
    align_ok = "PASS" in str(results["step3_minimap2"].get("status",""))
    results["step4_sv_calling"]  = simulate_sv_calling(resources)
    results["step5_medaka"]      = simulate_medaka(resources)
    results["step6_methylation"] = simulate_methylation(resources)

    weights = {"step1_basecalling":0.10,"step2_qc":0.05,"step3_minimap2":0.30,
               "step4_sv_calling":0.20,"step5_medaka":0.20,"step6_methylation":0.15}
    step_scores = {}
    for key in weights:
        s = results[key]
        if "FAIL" in str(s.get("status","")):
            step_scores[key] = 0
        else:
            step_scores[key] = min(100, s.get("throughput_reads_per_sec",0)/1000)

    weighted = sum(step_scores[k]*w for k,w in weights.items())
    fails    = sum(1 for s in results.values() if "FAIL" in str(s.get("status","")))

    if ram_gb < 16:
        score = round(min(weighted, 15), 1)
        capability = "🔴 Cannot run long read pipeline — RAM insufficient"
    elif fails >= 2:
        score = round(min(weighted, 35), 1)
        capability = "🟡 Partial pipeline only"
    else:
        score = round(weighted, 1)
        capability = "🟢 Full long read pipeline capable"

    total_time = round(sum(v.get("time_seconds",0) for v in results.values()), 2)

    print(f"\n  {'Step':<40} {'Status':<25} {'Time':>8}")
    print("  " + "-"*75)
    for key, r in results.items():
        name   = key.replace("_"," ").title()
        status = r.get("status","N/A")
        t      = r.get("time_seconds",0)
        print(f"  {name:<40} {status:<25} {t:>7.3f}s")
    print("  " + "-"*75)
    print(f"\n  ✅ Long Read Benchmark Complete!")
    print(f"  🏆 Score: {score}/100 | 💡 {capability}")

    return {"module":"LongRead","score":score,"capability":capability,
            "total_time_seconds":total_time,"hardware":resources,
            "hardware_warnings":warnings,"pipeline_steps":results}
