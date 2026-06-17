# BioMark - Genome Assembly Benchmark Module
# Simulates SPAdes, Flye, Hifiasm workloads
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
            "message":f"{ram}GB — Cannot run genome assemblers (need 16-500GB+)",
            "recommendation":"SPAdes needs 16GB min, Flye needs 32GB+, large genomes need 500GB+","can_run":False})
    elif ram < 32:
        warnings_list.append({"level":"🔴 ERROR","component":"RAM",
            "message":f"{ram}GB — Only tiny genomes (<50MB) possible",
            "recommendation":"32GB minimum for bacterial genomes, 64GB+ for human","can_run":False})
    elif ram < 64:
        warnings_list.append({"level":"🟡 WARNING","component":"RAM",
            "message":f"{ram}GB — Bacterial/small eukaryote genomes only",
            "recommendation":"64GB+ for human genome assembly","can_run":True})
    else:
        warnings_list.append({"level":"🟢 GOOD","component":"RAM",
            "message":f"{ram}GB — Suitable for most genome assemblies",
            "recommendation":"256GB+ for large genomes (>3GB)","can_run":True})

    if ssd < 200:
        warnings_list.append({"level":"🔴 CRITICAL","component":"SSD",
            "message":f"{ssd}GB free — Cannot store assembly intermediate files",
            "recommendation":"Assembly generates huge temp files — need 500GB+ free","can_run":False})
    else:
        warnings_list.append({"level":"🟢 GOOD","component":"SSD",
            "message":f"{ssd}GB free — Check space carefully","recommendation":"500GB+ recommended","can_run":True})

    ram_ok = ram >= 16
    warnings_list.append({
        "level": "⚪ N/A" if not ram_ok else "🟢 GOOD",
        "component": "CPU",
        "message": f"{resources['cpu_cores']} cores — {'Irrelevant while RAM insufficient' if not ram_ok else 'Good for assembly'}",
        "recommendation": "Fix RAM first" if not ram_ok else "16+ cores for faster assembly",
        "can_run": ram_ok
    })
    return warnings_list

def simulate_spades(resources):
    print("    🔧 Step 1: Short Read Assembly (SPAdes)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]
    ram_ok = ram_gb >= 16

    if not ram_ok:
        return {"step":"SPAdes","status":"🔴 FAIL",
                "note":f"SPAdes needs 16GB+ RAM. You have {ram_gb}GB.",
                "time_seconds":round(time.time()-start,3),"throughput_reads_per_sec":0}

    n_contigs  = random.randint(100, 5000)
    n50        = random.randint(10000, 5000000)
    total_len  = random.randint(1000000, 50000000)
    elapsed    = max(round(time.time()-start,3), 0.001)

    return {"step":"SPAdes","status":"✅ PASS","n_contigs":n_contigs,
            "n50_bp":n50,"total_assembly_bp":total_len,
            "largest_contig_bp":random.randint(n50,n50*5),
            "real_time_estimate":"1–24 hours depending on genome size",
            "time_seconds":elapsed,"throughput_reads_per_sec":safe_tp(n_contigs,elapsed)}

def simulate_flye(resources):
    print("    🔧 Step 2: Long Read Assembly (Flye)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]
    ram_ok = ram_gb >= 32

    if not ram_ok:
        return {"step":"Flye","status":"🔴 FAIL",
                "note":f"Flye needs 32GB+ RAM for human genome. You have {ram_gb}GB.",
                "time_seconds":round(time.time()-start,3),"throughput_reads_per_sec":0}

    n_contigs = random.randint(20, 500)
    n50       = random.randint(1000000, 50000000)
    elapsed   = max(round(time.time()-start,3), 0.001)

    return {"step":"Flye","status":"✅ PASS","n_contigs":n_contigs,
            "n50_bp":n50,"input":"Long reads (Nanopore/PacBio)",
            "real_time_estimate":"4–72 hours depending on genome size + coverage",
            "time_seconds":elapsed,"throughput_reads_per_sec":safe_tp(n_contigs,elapsed)}

def simulate_hifiasm(resources):
    print("    🔧 Step 3: HiFi Assembly (Hifiasm)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]
    ram_ok = ram_gb >= 64

    if not ram_ok:
        return {"step":"Hifiasm","status":"🔴 FAIL",
                "note":f"Hifiasm needs 64GB+ RAM for human genome. You have {ram_gb}GB.",
                "time_seconds":round(time.time()-start,3),"throughput_reads_per_sec":0}

    n_contigs = random.randint(40, 100)
    n50       = random.randint(10000000, 100000000)
    elapsed   = max(round(time.time()-start,3), 0.001)

    return {"step":"Hifiasm","status":"✅ PASS","n_contigs":n_contigs,
            "n50_bp":n50,"input":"PacBio HiFi reads","assembly_type":"Phased diploid",
            "real_time_estimate":"6–48 hours",
            "time_seconds":elapsed,"throughput_reads_per_sec":safe_tp(n_contigs,elapsed)}

def simulate_quast(resources):
    print("    📊 Step 4: Assembly QC (QUAST)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]
    if ram_gb < 16:
        return {"step":"QUAST","status":"🔴 FAIL","note":"No assembly to evaluate",
                "time_seconds":0,"throughput_reads_per_sec":0}

    elapsed = max(round(time.time()-start,3), 0.001)
    return {"step":"QUAST","status":"✅ PASS",
            "metrics":["N50","L50","Total length","# contigs","Largest contig"],
            "reference_comparison":"hg38","time_seconds":elapsed,
            "throughput_reads_per_sec":safe_tp(100,elapsed)}

def simulate_busco(resources):
    print("    🔍 Step 5: Completeness Check (BUSCO)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]
    if ram_gb < 16:
        return {"step":"BUSCO","status":"🔴 FAIL","note":"No assembly to evaluate",
                "time_seconds":0,"throughput_reads_per_sec":0}

    complete    = round(random.uniform(85, 99), 1)
    fragmented  = round(random.uniform(0.5, 5), 1)
    missing     = round(100 - complete - fragmented, 1)
    elapsed     = max(round(time.time()-start,3), 0.001)

    return {"step":"BUSCO","status":"✅ PASS",
            "complete_pct":complete,"fragmented_pct":fragmented,"missing_pct":missing,
            "lineage":"mammalia_odb10","time_seconds":elapsed,
            "throughput_reads_per_sec":safe_tp(1000,elapsed)}

def run_assembly_benchmark():
    print("\n🧬 BioMark — Genome Assembly Benchmark")
    print("   Tools : SPAdes + Flye + Hifiasm + QUAST + BUSCO")
    print("   Note  : Most RAM-intensive workflow in bioinformatics")
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
    results["step1_spades"]   = simulate_spades(resources)
    results["step2_flye"]     = simulate_flye(resources)
    results["step3_hifiasm"]  = simulate_hifiasm(resources)
    results["step4_quast"]    = simulate_quast(resources)
    results["step5_busco"]    = simulate_busco(resources)

    weights = {"step1_spades":0.25,"step2_flye":0.25,"step3_hifiasm":0.25,
               "step4_quast":0.15,"step5_busco":0.10}
    step_scores = {}
    for key in weights:
        s = results[key]
        if "FAIL" in str(s.get("status","")):
            step_scores[key] = 0
        else:
            step_scores[key] = min(100, s.get("throughput_reads_per_sec",0)/100)

    weighted = sum(step_scores[k]*w for k,w in weights.items())
    fails    = sum(1 for s in results.values() if "FAIL" in str(s.get("status","")))

    if ram_gb < 16:
        score = round(min(weighted, 5), 1)
        capability = "🔴 Cannot run any assembler — RAM critically insufficient"
    elif ram_gb < 32:
        score = round(min(weighted, 15), 1)
        capability = "🔴 Only tiny genomes possible — human genome requires 64GB+"
    elif ram_gb < 64:
        score = round(min(weighted, 40), 1)
        capability = "🟡 Bacterial genomes only — human genome needs 64GB+"
    else:
        score = round(weighted, 1)
        capability = "🟢 Full genome assembly capable"

    total_time = round(sum(v.get("time_seconds",0) for v in results.values()), 2)

    print(f"\n  {'Step':<40} {'Status':<25} {'Time':>8}")
    print("  " + "-"*75)
    for key, r in results.items():
        name   = key.replace("_"," ").title()
        status = r.get("status","N/A")
        t      = r.get("time_seconds",0)
        print(f"  {name:<40} {status:<25} {t:>7.3f}s")
    print("  " + "-"*75)
    print(f"\n  ✅ Assembly Benchmark Complete!")
    print(f"  🏆 Score: {score}/100 | 💡 {capability}")

    return {"module":"Assembly","score":score,"capability":capability,
            "total_time_seconds":total_time,"hardware":resources,
            "hardware_warnings":warnings,"pipeline_steps":results}
