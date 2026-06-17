# BioMark - Metagenomics Benchmark Module
# Simulates Kraken2, MetaPhlAn4, QIIME2 workloads
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

    if ram < 8:
        warnings_list.append({"level":"🔴 CRITICAL","component":"RAM",
            "message":f"{ram}GB — Cannot run Kraken2 (needs 8GB+ for database)",
            "recommendation":"Minimum 8GB RAM for Kraken2 standard database","can_run":False})
    elif ram < 16:
        warnings_list.append({"level":"🔴 ERROR","component":"RAM",
            "message":f"{ram}GB — Kraken2 standard DB needs 16GB RAM",
            "recommendation":"Use Kraken2 mini database (8GB) as workaround","can_run":False})
    elif ram < 32:
        warnings_list.append({"level":"🟡 WARNING","component":"RAM",
            "message":f"{ram}GB — Standard Kraken2 DB tight, use mini DB",
            "recommendation":"32GB for full Kraken2 standard database","can_run":True})
    else:
        warnings_list.append({"level":"🟢 GOOD","component":"RAM",
            "message":f"{ram}GB — Suitable for full metagenomics pipeline",
            "recommendation":"Sufficient for Kraken2 + MetaPhlAn4 + QIIME2","can_run":True})

    if ssd < 200:
        warnings_list.append({"level":"🔴 CRITICAL","component":"SSD",
            "message":f"{ssd}GB free — Cannot store Kraken2 database (~100GB)",
            "recommendation":"Need 200GB+ free. Kraken2 DB ~100GB, can use external SSD","can_run":False})
    else:
        warnings_list.append({"level":"🟢 GOOD","component":"SSD",
            "message":f"{ssd}GB free — Sufficient for metagenomics","recommendation":"Sufficient","can_run":True})

    ram_ok = ram >= 16
    warnings_list.append({
        "level": "⚪ N/A" if not ram_ok else "🟢 GOOD",
        "component": "CPU",
        "message": f"{resources['cpu_cores']} cores — {'Irrelevant while RAM insufficient' if not ram_ok else 'Good for metagenomics'}",
        "recommendation": "Fix RAM first" if not ram_ok else "8+ cores for faster classification",
        "can_run": ram_ok
    })
    return warnings_list

def simulate_qc(n_reads=100000):
    print("    📊 Step 1: QC (FastQC + fastp)...")
    start = time.time()
    passed = sum(1 for _ in range(n_reads) if random.uniform(0,40) >= 20)
    elapsed = max(round(time.time()-start,3), 0.001)
    return {"step":"QC","status":"✅ PASS","reads_passed":passed,
            "time_seconds":elapsed,"throughput_reads_per_sec":safe_tp(n_reads,elapsed)}

def simulate_kraken2(resources, n_reads=100000):
    print("    🦠 Step 2: Taxonomic Classification (Kraken2)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]
    ram_ok = ram_gb >= 16

    if not ram_ok:
        return {"step":"Kraken2","status":"🔴 FAIL",
                "note":f"Kraken2 standard DB needs 16GB+ RAM. You have {ram_gb}GB. Tip: Use mini database (8GB) as workaround.",
                "time_seconds":round(time.time()-start,3),"throughput_reads_per_sec":0}

    taxa = ["Bacteria","Archaea","Viruses","Eukaryota","Unclassified"]
    classified = int(n_reads * random.uniform(0.7, 0.95))
    results = []
    for i in range(50):
        results.append({
            "taxon"    : f"Species_{i:03d}",
            "kingdom"  : random.choice(taxa),
            "reads"    : random.randint(100, 10000),
            "confidence": round(random.uniform(0.5, 1.0), 3)
        })
    elapsed = max(round(time.time()-start,3), 0.001)
    return {
        "step"                   : "Kraken2",
        "status"                 : "✅ PASS",
        "reads_classified"       : classified,
        "classification_rate"    : round(classified/n_reads*100, 2),
        "taxa_identified"        : len(results),
        "database"               : "Kraken2 standard (~100GB)",
        "note"                   : "Database can be stored on external Thunderbolt SSD",
        "time_seconds"           : elapsed,
        "throughput_reads_per_sec": safe_tp(n_reads, elapsed)
    }

def simulate_bracken(resources, n_taxa=50):
    print("    📊 Step 3: Abundance Estimation (Bracken)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]
    if ram_gb < 16:
        return {"step":"Bracken","status":"🔴 FAIL",
                "note":"Depends on Kraken2 — skipped",
                "time_seconds":0,"throughput_reads_per_sec":0}
    abundances = [{"taxon":f"Species_{i}","abundance":round(random.uniform(0,1),4)} for i in range(n_taxa)]
    total = sum(a["abundance"] for a in abundances)
    for a in abundances:
        a["relative_abundance"] = round(a["abundance"]/total, 4)
    elapsed = max(round(time.time()-start,3), 0.001)
    return {"step":"Bracken","status":"✅ PASS","taxa_quantified":n_taxa,
            "dominant_taxon":abundances[0]["taxon"],
            "time_seconds":elapsed,"throughput_reads_per_sec":safe_tp(n_taxa,elapsed)}

def simulate_metaphlan(resources, n_reads=100000):
    print("    🧬 Step 4: Microbial Profiling (MetaPhlAn4)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]
    if ram_gb < 8:
        return {"step":"MetaPhlAn4","status":"🔴 FAIL",
                "note":f"MetaPhlAn4 needs 8GB+ RAM. You have {ram_gb}GB.",
                "time_seconds":round(time.time()-start,3),"throughput_reads_per_sec":0}
    species = [{"species":f"s__Species_{i}","relative_abundance":round(random.uniform(0,1),4)} for i in range(30)]
    elapsed = max(round(time.time()-start,3), 0.001)
    return {"step":"MetaPhlAn4","status":"✅ PASS","species_identified":len(species),
            "database":"ChocoPhlAn (~3GB)","time_seconds":elapsed,
            "throughput_reads_per_sec":safe_tp(n_reads,elapsed)}

def simulate_qiime2(resources, n_samples=20):
    print("    🌿 Step 5: 16S rRNA Analysis (QIIME2)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]
    if ram_gb < 8:
        return {"step":"QIIME2","status":"🔴 FAIL",
                "note":f"QIIME2 needs 8GB+ RAM. You have {ram_gb}GB.",
                "time_seconds":round(time.time()-start,3),"throughput_reads_per_sec":0}
    n_asvs         = random.randint(500, 3000)
    alpha_diversity = round(random.uniform(2, 6), 3)
    elapsed         = max(round(time.time()-start,3), 0.001)
    return {"step":"QIIME2","status":"✅ PASS","samples_processed":n_samples,
            "asvs_identified":n_asvs,"alpha_diversity_shannon":alpha_diversity,
            "time_seconds":elapsed,"throughput_reads_per_sec":safe_tp(n_samples,elapsed)}

def simulate_humann(resources, n_reads=100000):
    print("    🔬 Step 6: Functional Profiling (HUMAnN3)...")
    start  = time.time()
    ram_gb = resources["total_ram_gb"]
    if ram_gb < 16:
        return {"step":"HUMAnN3","status":"🔴 FAIL",
                "note":f"HUMAnN3 needs 16GB+ RAM. You have {ram_gb}GB.",
                "time_seconds":round(time.time()-start,3),"throughput_reads_per_sec":0}
    pathways = random.randint(100, 500)
    elapsed  = max(round(time.time()-start,3), 0.001)
    return {"step":"HUMAnN3","status":"✅ PASS","pathways_identified":pathways,
            "databases":["UniRef90","MetaCyc"],"time_seconds":elapsed,
            "throughput_reads_per_sec":safe_tp(n_reads,elapsed)}

def run_metagenomics_benchmark():
    print("\n🦠 BioMark — Metagenomics Benchmark")
    print("   Tools : Kraken2 + Bracken + MetaPhlAn4 + QIIME2 + HUMAnN3")
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
    results["step1_qc"]       = simulate_qc()
    results["step2_kraken2"]  = simulate_kraken2(resources)
    k2_ok = "PASS" in str(results["step2_kraken2"].get("status",""))
    results["step3_bracken"]  = simulate_bracken(resources)
    results["step4_metaphlan"]= simulate_metaphlan(resources)
    results["step5_qiime2"]   = simulate_qiime2(resources)
    results["step6_humann"]   = simulate_humann(resources)

    weights = {"step1_qc":0.05,"step2_kraken2":0.30,"step3_bracken":0.15,
               "step4_metaphlan":0.20,"step5_qiime2":0.15,"step6_humann":0.15}
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
        score = round(min(weighted, 20), 1)
        capability = "🔴 Cannot run full metagenomics — RAM insufficient"
    elif fails >= 2:
        score = round(min(weighted, 40), 1)
        capability = "🟡 Partial pipeline only"
    else:
        score = round(weighted, 1)
        capability = "🟢 Full metagenomics pipeline capable"

    total_time = round(sum(v.get("time_seconds",0) for v in results.values()), 2)

    print(f"\n  {'Step':<40} {'Status':<25} {'Time':>8}")
    print("  " + "-"*75)
    for key, r in results.items():
        name   = key.replace("_"," ").title()
        status = r.get("status","N/A")
        t      = r.get("time_seconds",0)
        print(f"  {name:<40} {status:<25} {t:>7.3f}s")
    print("  " + "-"*75)
    print(f"\n  ✅ Metagenomics Benchmark Complete!")
    print(f"  🏆 Score: {score}/100 | 💡 {capability}")

    return {"module":"Metagenomics","score":score,"capability":capability,
            "total_time_seconds":total_time,"hardware":resources,
            "hardware_warnings":warnings,"pipeline_steps":results}
