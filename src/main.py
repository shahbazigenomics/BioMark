# BioMark - Bioinformatics Benchmark Tool
# Author: Amir Shahbazi
# GitHub: shahbazigenomics

import argparse
import json
import os
import sys
import subprocess
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from system_info import get_full_system_profile, print_system_profile
from report import generate_html_report
from modules.DNA.wes_benchmark import run_dna_benchmark
from modules.RNA.rnaseq_benchmark import run_rnaseq_benchmark
from modules.RNA.scrna_benchmark import run_scrna_benchmark
from modules.Protein.protein_benchmark import run_protein_benchmark
from modules.Epigenomics.epigenomics_benchmark import run_epigenomics_benchmark
from modules.Metagenomics.metagenomics_benchmark import run_metagenomics_benchmark
from modules.LongRead.longread_benchmark import run_longread_benchmark
from modules.Assembly.assembly_benchmark import run_assembly_benchmark
from utils import BioMarkEncoder

def print_banner():
    print("""
╔══════════════════════════════════════════╗
║           BioMark v0.5.0  🧬             ║
║  Bioinformatics Benchmark Tool           ║
║  github.com/shahbazigenomics/BioMark     ║
╚══════════════════════════════════════════╝
    """)

def save_results(all_results, system_profile):
    os.makedirs("results", exist_ok=True)
    timestamp   = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = f"results/biomark_{timestamp}.json"
    final = {
        "system_profile"   : system_profile,
        "benchmark_results": all_results
    }
    with open(output_file, "w") as f:
        json.dump(final, f, indent=4, cls=BioMarkEncoder)
    print(f"\n💾 JSON saved to: {output_file}")
    return output_file

def print_summary(all_results):
    print("\n" + "=" * 55)
    print("📊 BioMark Final Score Summary")
    print("=" * 55)
    total_score = 0
    count       = 0
    for module, data in all_results.items():
        score  = data.get("score", "N/A")
        time_v = data.get("total_time_seconds", "N/A")
        print(f"  {module:<25} Score: {score}/100  Time: {time_v}s")
        if isinstance(score, (int, float)):
            total_score += score
            count       += 1
    if count > 0:
        overall = round(total_score / count, 1)
        print("-" * 55)
        print(f"  {'⭐ OVERALL SCORE':<25} {overall}/100")
    print("=" * 55)

def create_share_file(all_results, system_profile):
    """Create a clean shareable JSON file for community database"""
    import platform

    os.makedirs("results", exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Build clean share payload
    share_data = {
        "biomark_version" : "0.5.0",
        "submitted_at"    : timestamp,
        "machine" : {
            "os"          : system_profile.get("os", {}).get("system", "Unknown"),
            "os_version"  : system_profile.get("os", {}).get("macos_version", "Unknown"),
            "chip"        : system_profile.get("cpu", {}).get("brand",
                            system_profile.get("cpu", {}).get("brand", "Unknown")),
            "arch"        : system_profile.get("cpu", {}).get("architecture", "Unknown"),
            "cpu_cores"   : system_profile.get("cpu", {}).get("physical_cores", 0),
            "ram_gb"      : system_profile.get("ram", {}).get("total_gb", 0),
            "ssd_total_gb": system_profile.get("storage", {}).get("total_gb", 0),
            "ssd_free_gb" : system_profile.get("storage", {}).get("free_gb", 0),
            "gpu"         : system_profile.get("gpu", {}).get("model", "Unknown"),
            "mac_model"   : system_profile.get("cpu", {}).get("mac_model", "Unknown"),
        },
        "scores" : {},
        "capabilities": {}
    }

    # Add scores and capabilities
    for module, data in all_results.items():
        share_data["scores"][module]       = data.get("score", 0)
        share_data["capabilities"][module] = data.get("capability", "")

    # Calculate overall
    scores = [v for v in share_data["scores"].values() if isinstance(v, (int, float))]
    share_data["overall_score"] = round(sum(scores)/len(scores), 1) if scores else 0

    # Save share file
    share_file = f"results/biomark_share_{timestamp}.json"
    with open(share_file, "w") as f:
        json.dump(share_data, f, indent=4, cls=BioMarkEncoder)

    print(f"""
╔══════════════════════════════════════════════════════════╗
║           🌍 Share Your BioMark Score!                   ║
╚══════════════════════════════════════════════════════════╝

Your shareable score file has been saved to:
  {share_file}

To contribute to the BioMark community database:

  1. Open a GitHub Issue:
     https://github.com/shahbazigenomics/BioMark/issues/new

  2. Title: "Score Submission — [Your Machine] [RAM]GB"
     Example: "Score Submission — MacBook Pro M5 Pro 36GB"

  3. Attach your share file: {share_file}

Your scores:
  Overall: {share_data['overall_score']}/100
""")

    for module, score in share_data["scores"].items():
        print(f"  {module:<20} {score}/100")

    print(f"""
Thank you for contributing to BioMark! 🧬
Your data helps researchers choose the right hardware.
══════════════════════════════════════════════════════════
""")
    return share_file


def main():
    parser = argparse.ArgumentParser(
        description="BioMark — Bioinformatics Benchmark Tool"
    )
    parser.add_argument("--all",          action="store_true", help="Run all benchmarks")
    parser.add_argument("--dna",          action="store_true", help="Run DNA/WES benchmark")
    parser.add_argument("--rna",          action="store_true", help="Run bulk RNA-seq benchmark")
    parser.add_argument("--scrna",        action="store_true", help="Run scRNA-seq benchmark")
    parser.add_argument("--protein",      action="store_true", help="Run Protein structure benchmark")
    parser.add_argument("--epigenomics",  action="store_true", help="Run Epigenomics benchmark")
    parser.add_argument("--metagenomics", action="store_true", help="Run Metagenomics benchmark")
    parser.add_argument("--longread",     action="store_true", help="Run Long Read benchmark")
    parser.add_argument("--assembly",     action="store_true", help="Run Genome Assembly benchmark")
    parser.add_argument("--sysinfo",      action="store_true", help="Show system hardware info only")
    parser.add_argument("--share",        action="store_true", help="Generate shareable score file for community database")
    args = parser.parse_args()

    print_banner()

    system_profile = get_full_system_profile()
    print_system_profile(system_profile)

    if args.sysinfo:
        return

    all_results = {}

    if args.dna or args.all:
        all_results["DNA_WES"]      = run_dna_benchmark()
    if args.rna or args.all:
        all_results["RNA_seq"]      = run_rnaseq_benchmark()
    if args.scrna or args.all:
        all_results["scRNA_seq"]    = run_scrna_benchmark()
    if args.protein or args.all:
        all_results["Protein"]      = run_protein_benchmark()
    if args.epigenomics or args.all:
        all_results["Epigenomics"]  = run_epigenomics_benchmark()
    if args.metagenomics or args.all:
        all_results["Metagenomics"] = run_metagenomics_benchmark()
    if args.longread or args.all:
        all_results["LongRead"]     = run_longread_benchmark()
    if args.assembly or args.all:
        all_results["Assembly"]     = run_assembly_benchmark()

    if not all_results:
        print("⚠️  Please specify a module:")
        print("  python src/main.py --dna")
        print("  python src/main.py --rna")
        print("  python src/main.py --scrna")
        print("  python src/main.py --protein")
        print("  python src/main.py --epigenomics")
        print("  python src/main.py --metagenomics")
        print("  python src/main.py --longread")
        print("  python src/main.py --assembly")
        print("  python src/main.py --all")
        print("  python src/main.py --sysinfo")
        return

    print_summary(all_results)
    save_results(all_results, system_profile)

    report_file = generate_html_report(system_profile, all_results)
    subprocess.run(["open", report_file])

    if args.share:
        create_share_file(all_results, system_profile)

    print("\n✅ BioMark benchmark complete!\n")

if __name__ == "__main__":
    main()
