# BioMark - WES Pipeline Benchmark Module
# Simulates real WES computational workloads
# Based on standard WES pipeline for rare disease analysis
# Author: Amir Shahbazi
# GitHub: shahbazigenomics

import time
import random
import os
import gzip
import multiprocessing
import psutil
import numpy as np

# ════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ════════════════════════════════════════════════════════════

def generate_dna_sequence(length):
    return ''.join(random.choices(['A','T','G','C'], k=length))

def generate_quality_scores(length):
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

def evaluate_hardware(resources):
    warnings = []
    ram   = resources["total_ram_gb"]
    ssd   = resources["free_ssd_gb"]
    cores = resources["cpu_cores"]

    # RAM
    if ram < 8:
        warnings.append({
            "level": "🔴 CRITICAL", "component": "RAM",
            "message": f"{ram}GB — Cannot run any WES analysis",
            "recommendation": "Minimum 16GB RAM required",
            "can_run": False
        })
    elif ram < 16:
        warnings.append({
            "level": "🔴 CRITICAL", "component": "RAM",
            "message": f"{ram}GB — Cannot load hg38 index (BWA-MEM2 needs 16GB+)",
            "recommendation": "Need 16GB+ for alignment step",
            "can_run": False
        })
    elif ram < 24:
        warnings.append({
            "level": "🟡 WARNING", "component": "RAM",
            "message": f"{ram}GB — Single sample only, alignment will be slow",
            "recommendation": "36GB recommended for comfortable WES",
            "can_run": True
        })
    elif ram < 36:
        warnings.append({
            "level": "🟡 CAUTION", "component": "RAM",
            "message": f"{ram}GB — Single sample ok, trio analysis tight",
            "recommendation": "36GB+ for trio analysis",
            "can_run": True
        })
    else:
        warnings.append({
            "level": "🟢 GOOD", "component": "RAM",
            "message": f"{ram}GB — Suitable for WES + trio analysis",
            "recommendation": "64GB+ for large cohort analysis",
            "can_run": True
        })

    # SSD
    if ssd < 100:
        warnings.append({
            "level": "🔴 CRITICAL", "component": "SSD",
            "message": f"{ssd}GB free — Cannot store reference + sample data",
            "recommendation": "Need 500GB+ free. Tip: gnomAD (~1TB) can be stored on external SSD (USB-C/Thunderbolt recommended for speed)",
            "can_run": False
        })
    elif ssd < 500:
        warnings.append({
            "level": "🔴 ERROR", "component": "SSD",
            "message": f"{ssd}GB free — Very tight for WES pipeline",
            "recommendation": "hg38(3GB)+FASTQ(10GB)+BAM(30GB) needed locally. gnomAD(~1TB) can be on external SSD via Thunderbolt",
            "can_run": False
        })
    elif ssd < 1000:
        warnings.append({
            "level": "🟡 WARNING", "component": "SSD",
            "message": f"{ssd}GB free — One sample at a time only",
            "recommendation": "1TB+ internal recommended. gnomAD(~1TB) can be stored on fast external SSD (Thunderbolt 3/4/5)",
            "can_run": True
        })
    else:
        warnings.append({
            "level": "🟢 GOOD", "component": "SSD",
            "message": f"{ssd}GB free — Sufficient for WES pipeline",
            "recommendation": "For gnomAD(~1TB) database consider a fast external Thunderbolt SSD if internal space is limited",
            "can_run": True
        })

    # CPU — only meaningful if RAM is sufficient
    ram_ok = ram >= 16
    if not ram_ok:
        warnings.append({
            "level": "⚪ N/A", "component": "CPU",
            "message": f"{cores} cores — CPU is irrelevant while RAM is insufficient",
            "recommendation": "Fix RAM first — CPU cannot help if hg38 index cannot load",
            "can_run": False
        })
    elif cores < 8:
        warnings.append({
            "level": "🟡 WARNING", "component": "CPU",
            "message": f"{cores} cores — Slow but functional",
            "recommendation": "8+ cores recommended",
            "can_run": True
        })
    elif cores < 16:
        warnings.append({
            "level": "🟢 GOOD", "component": "CPU",
            "message": f"{cores} cores — Good for single sample WES",
            "recommendation": "16+ cores for faster trio analysis",
            "can_run": True
        })
    else:
        warnings.append({
            "level": "🟢 GREAT", "component": "CPU",
            "message": f"{cores} cores — Excellent for WES + trio",
            "recommendation": "Ideal for parallel multi-sample processing",
            "can_run": True
        })

    return warnings


# ════════════════════════════════════════════════════════════
# PIPELINE STEPS
# ════════════════════════════════════════════════════════════

def simulate_fastqc(n_reads=100000, read_length=150):
    print("    📊 Step 1: QC (FastQC)...")
    start = time.time()
    os.makedirs("results", exist_ok=True)
    test_file = "results/biomark_qc_test.fastq.gz"

    with gzip.open(test_file, 'wt') as f:
        for i in range(n_reads):
            seq  = generate_dna_sequence(read_length)
            qual = generate_quality_scores(read_length)
            f.write(f"@READ_{i}\n{seq}\n+\n{qual}\n")

    file_size_mb = os.path.getsize(test_file) / (1024*1024)
    gc_list = []
    reads_parsed = 0
    read_start = time.time()

    with gzip.open(test_file, 'rt') as f:
        while True:
            if not f.readline(): break
            seq  = f.readline().strip()
            f.readline(); f.readline()
            gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
            gc_list.append(gc)
            reads_parsed += 1

    read_time = round(time.time() - read_start, 3)
    os.remove(test_file)
    elapsed = round(time.time() - start, 3)

    mean_gc   = round(sum(gc_list)/len(gc_list), 2) if gc_list else 0
    speed     = round(file_size_mb / read_time, 2)
    throughput = round(reads_parsed / read_time)

    # Step scoring
    # QC only needs SSD speed — any machine can do this
    step_score = min(100, speed / 5)  # 500MB/s = 100 score
    status = "✅ PASS"

    return {
        "step"          : "QC_FastQC",
        "status"        : status,
        "step_score"    : round(step_score, 1),
        "reads_analyzed": reads_parsed,
        "mean_gc"       : mean_gc,
        "read_speed_mbs": speed,
        "bottleneck"    : "SSD read speed",
        "time_seconds"  : elapsed,
        "throughput_reads_per_sec": throughput
    }


def simulate_trimmomatic(n_reads=100000, read_length=150):
    print("    ✂️  Step 2: Trimming (Trimmomatic)...")
    start = time.time()
    passed = trimmed = dropped = 0

    for i in range(n_reads):
        seq  = generate_dna_sequence(read_length)
        qual = generate_quality_scores(read_length)
        tlen = read_length
        for j in range(0, len(qual)-4, 4):
            wq = sum(ord(c)-33 for c in qual[j:j+4]) / 4
            if wq < 15:
                tlen = j
                break
        if tlen >= 36:
            passed += 1
            if tlen < read_length:
                trimmed += 1
        else:
            dropped += 1

    elapsed    = round(time.time() - start, 3)
    throughput = round(n_reads / elapsed)
    step_score = min(100, throughput / 50000)
    status     = "✅ PASS"

    return {
        "step"          : "Trimming_Trimmomatic",
        "status"        : status,
        "step_score"    : round(step_score, 1),
        "reads_passed"  : passed,
        "reads_dropped" : dropped,
        "survival_rate" : round(passed/n_reads*100, 2),
        "bottleneck"    : "CPU speed",
        "time_seconds"  : elapsed,
        "throughput_reads_per_sec": throughput
    }


def simulate_alignment_bwamem2(n_reads=50000, read_length=150):
    print("    🔗 Step 3a: Alignment (BWA-MEM2)...")
    start     = time.time()
    resources = get_system_resources()
    ram_gb    = resources["total_ram_gb"]

    # CRITICAL CHECK: BWA-MEM2 needs 16GB RAM for hg38 index
    ram_sufficient = ram_gb >= 16
    if not ram_sufficient:
        status     = "🔴 FAIL — Insufficient RAM"
        step_score = 0
        note       = f"BWA-MEM2 needs 16GB+ RAM for hg38 index. You have {ram_gb}GB."
        elapsed    = round(time.time() - start, 3)
        return {
            "step"          : "Alignment_BWA_MEM2",
            "status"        : status,
            "step_score"    : step_score,
            "note"          : note,
            "ram_required"  : "16GB",
            "ram_available" : f"{ram_gb}GB",
            "bottleneck"    : "RAM (critical)",
            "time_seconds"  : elapsed,
            "throughput_reads_per_sec": 0
        }

    # If RAM ok, simulate alignment
    try:
        alloc = min(
            300_000_000,
            int(resources["available_ram_gb"] * 0.3 * 1024**3 / 8)
        )
        ref = np.random.randint(0, 4, size=alloc, dtype=np.uint8)
    except MemoryError:
        ref = None

    aligned = multi = unmapped = 0
    align_start = time.time()

    for i in range(n_reads):
        read  = generate_dna_sequence(read_length)
        score = random.randint(90, 150) if ref is not None else random.randint(60, 150)
        if score >= 120:
            aligned += 1
        elif score >= 90:
            multi += 1
        else:
            unmapped += 1

    if ref is not None:
        del ref

    align_time = round(time.time() - align_start, 3)
    elapsed    = round(time.time() - start, 3)
    throughput = round(n_reads / align_time) if align_time > 0 else 0
    step_score = min(100, throughput / 5000)
    status     = "✅ PASS"

    return {
        "step"          : "Alignment_BWA_MEM2",
        "status"        : status,
        "step_score"    : round(step_score, 1),
        "reads_aligned" : aligned,
        "unmapped"      : unmapped,
        "multi_mapped"  : multi,
        "mapping_rate"  : round(aligned/n_reads*100, 2),
        "ram_required"  : "16GB",
        "ram_available" : f"{ram_gb}GB",
        "bottleneck"    : "RAM bandwidth + CPU cores",
        "time_seconds"  : elapsed,
        "throughput_reads_per_sec": throughput
    }


def simulate_alignment_hisat2(n_reads=50000, read_length=150):
    print("    🔗 Step 3b: Alignment (HISAT2)...")
    start     = time.time()
    resources = get_system_resources()
    ram_gb    = resources["total_ram_gb"]

    # HISAT2 needs ~8GB for hg38
    ram_sufficient = ram_gb >= 8
    if not ram_sufficient:
        return {
            "step"      : "Alignment_HISAT2",
            "status"    : "🔴 FAIL — Insufficient RAM",
            "step_score": 0,
            "note"      : f"HISAT2 needs 8GB+ RAM. You have {ram_gb}GB.",
            "time_seconds": round(time.time()-start, 3),
            "throughput_reads_per_sec": 0
        }

    aligned = multi = unmapped = 0
    for i in range(n_reads):
        score = random.gauss(128, 18)
        if score >= 120:
            aligned += 1
        elif score >= 85:
            multi += 1
        else:
            unmapped += 1

    elapsed    = round(time.time() - start, 3)
    throughput = round(n_reads / elapsed)
    step_score = min(100, throughput / 3000)

    return {
        "step"          : "Alignment_HISAT2",
        "status"        : "✅ PASS",
        "step_score"    : round(step_score, 1),
        "reads_aligned" : aligned,
        "mapping_rate"  : round(aligned/n_reads*100, 2),
        "ram_required"  : "8GB",
        "ram_available" : f"{ram_gb}GB",
        "bottleneck"    : "CPU cores",
        "time_seconds"  : elapsed,
        "throughput_reads_per_sec": throughput
    }


def simulate_post_alignment(n_reads=100000):
    print("    🔧 Step 4: Post-alignment (Samtools + GATK)...")
    start   = time.time()
    results = {}

    # 4a Filter unmapped
    s = time.time()
    unmapped = sum(
        1 for _ in range(n_reads) if random.randint(0,2048) & 4
    )
    results["filter_unmapped"] = {
        "removed": unmapped,
        "time_seconds": round(time.time()-s, 3)
    }

    # 4b Filter multi-mapped
    s = time.time()
    multi = sum(
        1 for _ in range(n_reads-unmapped)
        if random.randint(0,60) < 20
    )
    results["filter_multimapped"] = {
        "removed": multi,
        "time_seconds": round(time.time()-s, 3)
    }

    # 4c Add read groups
    s = time.time()
    remaining = n_reads - unmapped - multi
    results["add_read_groups"] = {
        "reads_processed": remaining,
        "time_seconds": round(time.time()-s, 3)
    }

    # 4d Mark duplicates
    s = time.time()
    positions = {}
    dups = 0
    for i in range(remaining):
        key = f"chr{random.randint(1,22)}:{random.randint(1,250_000_000)}"
        if key in positions:
            dups += 1
        else:
            positions[key] = i
    results["mark_duplicates"] = {
        "duplicates_marked": dups,
        "duplication_rate" : round(dups/remaining*100, 2) if remaining > 0 else 0,
        "time_seconds"     : round(time.time()-s, 3)
    }

    # 4e BQSR
    s = time.time()
    after_dedup = remaining - dups
    recal = 0
    for i in range(min(after_dedup, 20000)):
        _ = min(40, int(random.randint(2,40) * random.uniform(0.8,1.2)))
        recal += 1
    results["bqsr"] = {
        "reads_recalibrated": recal,
        "time_seconds": round(time.time()-s, 3)
    }

    elapsed    = round(time.time() - start, 3)
    throughput = round(n_reads / elapsed)
    step_score = min(100, throughput / 20000)

    return {
        "step"          : "Post_Alignment",
        "status"        : "✅ PASS",
        "step_score"    : round(step_score, 1),
        "substeps"      : results,
        "reads_input"   : n_reads,
        "reads_output"  : after_dedup,
        "bottleneck"    : "RAM + SSD I/O",
        "time_seconds"  : elapsed,
        "throughput_reads_per_sec": throughput
    }


def _call_variants(n_positions, tool):
    variants = 0
    chromosomes = [f"chr{i}" for i in range(1,23)] + ["chrX","chrY"]
    for i in range(n_positions):
        depth = random.randint(10, 200)
        alt   = random.randint(0, depth)
        vaf   = alt/depth if depth > 0 else 0
        qual  = random.uniform(0, 100)
        if tool == "GATK":
            if vaf >= 0.2 and qual >= 30 and depth >= 10:
                variants += 1
        elif tool == "bcftools":
            if vaf >= 0.15 and qual >= 20 and depth >= 8:
                variants += 1
        elif tool == "FreeBayes":
            if vaf >= 0.1 and qual >= 20 and depth >= 6:
                variants += 1
    return variants

def simulate_variant_calling_gatk(n_positions=200000):
    print("    🧬 Step 5a: Variant Calling (GATK)...")
    start    = time.time()
    variants = _call_variants(n_positions, "GATK")
    elapsed  = round(time.time() - start, 3)
    tp       = round(n_positions / elapsed)
    return {
        "step"          : "Variant_Calling_GATK",
        "status"        : "✅ PASS",
        "step_score"    : min(100, round(tp/200000, 1)),
        "tool"          : "GATK HaplotypeCaller",
        "variants_called": variants,
        "bottleneck"    : "CPU + RAM",
        "time_seconds"  : elapsed,
        "throughput_reads_per_sec": tp
    }

def simulate_variant_calling_bcftools(n_positions=200000):
    print("    🧬 Step 5b: Variant Calling (bcftools)...")
    start    = time.time()
    variants = _call_variants(n_positions, "bcftools")
    elapsed  = round(time.time() - start, 3)
    tp       = round(n_positions / elapsed)
    return {
        "step"          : "Variant_Calling_bcftools",
        "status"        : "✅ PASS",
        "step_score"    : min(100, round(tp/200000, 1)),
        "tool"          : "bcftools mpileup",
        "variants_called": variants,
        "bottleneck"    : "CPU + SSD I/O",
        "time_seconds"  : elapsed,
        "throughput_reads_per_sec": tp
    }

def simulate_variant_calling_freebayes(n_positions=200000):
    print("    🧬 Step 5c: Variant Calling (FreeBayes)...")
    start    = time.time()
    variants = _call_variants(n_positions, "FreeBayes")
    elapsed  = round(time.time() - start, 3)
    tp       = round(n_positions / elapsed)
    return {
        "step"          : "Variant_Calling_FreeBayes",
        "status"        : "✅ PASS",
        "step_score"    : min(100, round(tp/200000, 1)),
        "tool"          : "FreeBayes",
        "variants_called": variants,
        "bottleneck"    : "CPU",
        "time_seconds"  : elapsed,
        "throughput_reads_per_sec": tp
    }


def simulate_annotation(n_variants=50000):
    print("    📚 Step 6: Annotation (ANNOVAR)...")
    start     = time.time()
    resources = get_system_resources()
    ssd_free  = resources["free_ssd_gb"]

    # gnomAD needs ~1TB SSD
    ssd_ok = ssd_free >= 1000
    if not ssd_ok:
        status = "⚠️  PARTIAL — gnomAD database too large for available SSD"
        note   = f"gnomAD ~1TB needed. Only {ssd_free}GB free. Local annotation limited."
    else:
        status = "✅ PASS"
        note   = "Full annotation with gnomAD possible"

    consequences = [
        "missense_variant","synonymous_variant","stop_gained",
        "frameshift_variant","splice_region_variant","intron_variant"
    ]
    annotated = 0
    for i in range(n_variants):
        _ = {
            "gnomad_af"  : random.uniform(0, 0.5),
            "cadd_score" : random.uniform(0, 40),
            "clinvar"    : random.choice(["Pathogenic","Benign","VUS"]),
            "consequence": random.choice(consequences)
        }
        annotated += 1

    elapsed = round(time.time() - start, 3)
    tp      = round(n_variants / elapsed)

    return {
        "step"          : "Annotation_ANNOVAR",
        "status"        : status,
        "step_score"    : min(100, round(tp/50000, 1)) if ssd_ok else 30,
        "note"          : note,
        "databases"     : ["gnomAD (~1TB)","CADD","ClinVar","dbSNP","OMIM"],
        "variants_annotated": annotated,
        "ssd_sufficient_for_gnomad": ssd_ok,
        "bottleneck"    : "SSD random read speed (gnomAD)",
        "time_seconds"  : elapsed,
        "throughput_reads_per_sec": tp
    }


def simulate_filtering(n_variants=50000):
    print("    🔍 Step 7: Variant Filtering...")
    start  = time.time()
    total  = n_variants
    after_maf = after_cadd = after_qual = after_func = final = 0

    for i in range(n_variants):
        af   = random.uniform(0, 0.5)
        cadd = random.uniform(0, 40)
        qual = random.uniform(0, 100)
        dp   = random.randint(1, 500)
        cons = random.choice([
            "missense_variant","synonymous_variant",
            "stop_gained","frameshift_variant",
            "splice_region_variant","intron_variant"
        ])
        if af >= 0.01: continue
        after_maf += 1
        if cadd < 20: continue
        after_cadd += 1
        if qual < 30 or dp < 10: continue
        after_qual += 1
        if cons not in [
            "missense_variant","stop_gained",
            "frameshift_variant","splice_region_variant"
        ]: continue
        after_func += 1
        final += 1

    elapsed = round(time.time() - start, 3)
    tp      = round(n_variants / elapsed)

    return {
        "step"            : "Variant_Filtering",
        "status"          : "✅ PASS",
        "step_score"      : min(100, round(tp/50000, 1)),
        "filters"         : ["MAF<0.01 (gnomAD)","CADD>=20","QUAL>=30 DP>=10","Functional consequence"],
        "variants_input"  : total,
        "after_maf"       : after_maf,
        "after_cadd"      : after_cadd,
        "after_quality"   : after_qual,
        "after_consequence": after_func,
        "final_candidates": final,
        "reduction_rate"  : round((1-final/total)*100, 2),
        "bottleneck"      : "CPU + database lookup",
        "time_seconds"    : elapsed,
        "throughput_reads_per_sec": tp
    }


# ════════════════════════════════════════════════════════════
# MAIN WES BENCHMARK
# ════════════════════════════════════════════════════════════

def run_dna_benchmark(mode="single"):
    print("\n🧬 BioMark — WES Pipeline Benchmark")
    print(f"   Mode     : {'Single Sample' if mode=='single' else 'Trio (3 samples)'}")
    print("   Reference: hg38 (GRCh38)")
    print("   Input    : ~10GB FASTQ, 100x depth, 150bp PE")
    print("=" * 55)

    resources = get_system_resources()
    warnings  = evaluate_hardware(resources)
    ram_gb    = resources["total_ram_gb"]
    ssd_free  = resources["free_ssd_gb"]

    print("\n  📋 Hardware Assessment:")
    for w in warnings:
        print(f"     {w['level']} [{w['component']}] {w['message']}")
        print(f"     💡 {w['recommendation']}")
    print()

    results = {}
    results["step1_qc"]            = simulate_fastqc()
    results["step2_trimming"]      = simulate_trimmomatic()
    results["step3a_bwa"]          = simulate_alignment_bwamem2()
    results["step3b_hisat2"]       = simulate_alignment_hisat2()
    results["step4_post"]          = simulate_post_alignment()
    results["step5a_gatk"]         = simulate_variant_calling_gatk()
    results["step5b_bcftools"]     = simulate_variant_calling_bcftools()
    results["step5c_freebayes"]    = simulate_variant_calling_freebayes()
    results["step6_annotation"]    = simulate_annotation()
    results["step7_filtering"]     = simulate_filtering()

    # ── HONEST SCORING ──
    # Each step has a step_score (0-100)
    # BUT: if a critical step fails → overall score is capped

    step_scores = {k: v.get("step_score", 0) for k,v in results.items()}

    # Check critical failures
    bwa_failed  = results["step3a_bwa"]["step_score"] == 0
    ssd_failed  = ssd_free < 100

    # Weighted average of step scores
    weights = {
        "step1_qc"         : 0.05,
        "step2_trimming"   : 0.05,
        "step3a_bwa"       : 0.20,  # Most important
        "step3b_hisat2"    : 0.10,
        "step4_post"       : 0.10,
        "step5a_gatk"      : 0.20,  # Most important
        "step5b_bcftools"  : 0.05,
        "step5c_freebayes" : 0.05,
        "step6_annotation" : 0.10,
        "step7_filtering"  : 0.10,
    }

    weighted_score = sum(
        step_scores.get(k, 0) * w
        for k, w in weights.items()
    )

    # Cap score if critical steps fail
    if bwa_failed and ssd_failed:
        # Both alignment and storage fail → very low score
        dna_score = round(min(weighted_score, 15), 1)
        capability = "🔴 Cannot run real WES pipeline (RAM + SSD insufficient)"
    elif bwa_failed:
        # Alignment fails → cap at 30
        dna_score  = round(min(weighted_score, 30), 1)
        capability = "🔴 Cannot run alignment step (RAM insufficient for hg38)"
    elif ssd_failed:
        dna_score  = round(min(weighted_score, 35), 1)
        capability = "🔴 Cannot store WES data (SSD space insufficient)"
    elif ram_gb < 24:
        dna_score  = round(min(weighted_score, 65), 1)
        capability = "🟡 Single sample WES only — no parallel processing"
    elif ram_gb < 36:
        dna_score  = round(weighted_score, 1)
        capability = "🟡 Single sample comfortable — trio analysis limited"
    else:
        dna_score  = round(weighted_score, 1)
        capability = "🟢 Full WES pipeline — single + trio capable"

    total_time = round(sum(
        v.get("time_seconds", 0) for v in results.values()
    ), 2)

    # ── Print Summary Table ──
    step_names = {
        "step1_qc"         : "Step 1 : QC (FastQC)",
        "step2_trimming"   : "Step 2 : Trimming (Trimmomatic)",
        "step3a_bwa"       : "Step 3a: Alignment (BWA-MEM2)",
        "step3b_hisat2"    : "Step 3b: Alignment (HISAT2)",
        "step4_post"       : "Step 4 : Post-alignment (GATK)",
        "step5a_gatk"      : "Step 5a: Variant Calling (GATK)",
        "step5b_bcftools"  : "Step 5b: Variant Calling (bcftools)",
        "step5c_freebayes" : "Step 5c: Variant Calling (FreeBayes)",
        "step6_annotation" : "Step 6 : Annotation (ANNOVAR)",
        "step7_filtering"  : "Step 7 : Variant Filtering",
    }

    print(f"\n  {'Step':<38} {'Status':<30} {'Score':>6}  {'Time':>8}")
    print("  " + "-" * 85)
    for key, name in step_names.items():
        if key in results:
            r      = results[key]
            status = r.get("status", "N/A")
            score  = r.get("step_score", 0)
            t      = r.get("time_seconds", 0)
            print(f"  {name:<38} {status:<30} {score:>5.1f}  {t:>7.3f}s")

    print("  " + "-" * 85)
    print(f"\n  💡 Capability: {capability}")
    print(f"\n  ✅ WES Benchmark Complete!")
    print(f"  ⏱  Total time    : {total_time}s")
    print(f"  🧠 RAM           : {ram_gb}GB")
    print(f"  💾 SSD free      : {ssd_free}GB")
    print(f"  ⚡ CPU cores     : {resources['cpu_cores']}")
    print(f"  🏆 WES Score     : {dna_score}/100")

    return {
        "module"            : "DNA_WES",
        "mode"              : mode,
        "score"             : dna_score,
        "capability"        : capability,
        "total_time_seconds": total_time,
        "hardware"          : resources,
        "hardware_warnings" : warnings,
        "pipeline_steps"    : results
    }
