# BioMark - HTML Report Generator
# Author: Amir Shahbazi
# GitHub: shahbazigenomics

import json
import os
from datetime import datetime


def _score_color(score):
    if not isinstance(score, (int, float)):
        return "#6b7280"
    if score >= 80: return "#16a34a"
    if score >= 50: return "#d97706"
    return "#dc2626"

def _score_badge(score):
    if not isinstance(score, (int, float)):
        return "N/A"
    if score >= 80: return f"🟢 {score}/100"
    if score >= 50: return f"🟡 {score}/100"
    return f"🔴 {score}/100"

def _build_wes_table(module_data):
    """Build detailed WES pipeline table"""
    rows = ""
    steps = module_data.get("pipeline_steps", {})

    step_labels = {
        "step1_qc"         : ("📊", "Step 1", "QC", "FastQC"),
        "step2_trimming"   : ("✂️",  "Step 2", "Trimming", "Trimmomatic"),
        "step3a_bwa"       : ("🔗", "Step 3a", "Alignment", "BWA-MEM2"),
        "step3b_hisat2"    : ("🔗", "Step 3b", "Alignment", "HISAT2"),
        "step4_post"       : ("🔧", "Step 4", "Post-alignment", "Samtools+GATK"),
        "step5a_gatk"      : ("🧬", "Step 5a", "Variant Calling", "GATK HaplotypeCaller"),
        "step5b_bcftools"  : ("🧬", "Step 5b", "Variant Calling", "bcftools"),
        "step5c_freebayes" : ("🧬", "Step 5c", "Variant Calling", "FreeBayes"),
        "step6_annotation" : ("📚", "Step 6", "Annotation", "ANNOVAR"),
        "step7_filtering"  : ("🔍", "Step 7", "Variant Filtering", "Custom filters"),
    }

    for key, (icon, step_num, step_name, tool) in step_labels.items():
        if key not in steps:
            continue
        s          = steps[key]
        status     = s.get("status", "N/A")
        score      = s.get("step_score", 0)
        t          = s.get("time_seconds", 0)
        tp         = s.get("throughput_reads_per_sec", 0)
        bottleneck = s.get("bottleneck", "N/A")
        note       = s.get("note", "")

        # Status styling
        if "FAIL" in status or "CRITICAL" in status:
            status_html = f'<span class="badge-fail">{status}</span>'
            row_style   = "background:#fff5f5"
        elif "PARTIAL" in status or "WARNING" in status:
            status_html = f'<span class="badge-warn">{status}</span>'
            row_style   = "background:#fffbeb"
        else:
            status_html = f'<span class="badge-pass">{status}</span>'
            row_style   = ""

        tp_str    = f"{tp:,}/s" if tp else "N/A"
        note_html = f'<br><small style="color:#6b7280">{note}</small>' if note else ""

        rows += f"""
        <tr style="{row_style}">
            <td style="font-size:1.1rem;text-align:center">{icon}</td>
            <td style="text-align:center"><strong>{step_num}</strong></td>
            <td style="text-align:center">{step_name}</td>
            <td style="text-align:center"><code>{tool}</code></td>
            <td>{status_html}{note_html}</td>
            <td style="color:#475569">{bottleneck}</td>
            <td style="color:#475569">{tp_str}</td>
            <td style="color:#475569">{t:.3f}s</td>
        </tr>"""

    return rows


def _build_rna_table(module_data):
    """Build detailed RNA-seq pipeline table with PASS/FAIL status"""
    rows  = ""
    steps = module_data.get("pipeline_steps", {})

    step_labels = {
        "step1_qc"             : ("📊", "Step 1",  "QC",               "FastQC + MultiQC"),
        "step2_trimming"       : ("✂️",  "Step 2",  "Trimming",         "Trimmomatic / fastp"),
        "step3a_star"          : ("🔗", "Step 3a", "Alignment",        "STAR (splice-aware)"),
        "step3b_hisat2"        : ("🔗", "Step 3b", "Alignment",        "HISAT2 (splice-aware)"),
        "step4_postprocess"    : ("🔧", "Step 4",  "Post-processing",  "Samtools"),
        "step5a_featurecounts" : ("🔢", "Step 5a", "Quantification",   "featureCounts"),
        "step5b_salmon"        : ("🔢", "Step 5b", "Quantification",   "Salmon"),
        "step6a_deseq2"        : ("📈", "Step 6a", "Diff. Expression", "DESeq2"),
        "step6b_edger"         : ("📈", "Step 6b", "Diff. Expression", "edgeR"),
        "step7_pathway"        : ("🧪", "Step 7",  "Pathway Analysis", "clusterProfiler"),
    }

    for key, (icon, step_num, step_name, tool) in step_labels.items():
        if key not in steps:
            continue

        s          = steps[key]
        t          = s.get("time_seconds", 0)
        tp         = s.get("throughput_reads_per_sec", 0)
        tp_str     = f"{tp:,}/s" if tp else "N/A"
        step_status = str(s.get("status", ""))
        star_ok    = s.get("star_ram_sufficient", True)

        # ── Determine PASS/FAIL ──
        if key == "step3a_star" and not star_ok:
            badge    = '<span class="badge-fail">🔴 FAIL</span>'
            row_bg   = "background:#fff5f5"
            extra    = f'<br><small style="color:#dc2626">STAR needs 32GB RAM — you have {s.get("star_ram_available_gb","?")}GB</small>'
        elif "FAIL" in step_status:
            badge    = '<span class="badge-fail">🔴 FAIL</span>'
            row_bg   = "background:#fff5f5"
            note     = s.get("note","")
            extra    = f'<br><small style="color:#dc2626">{note}</small>' if note else ""
        elif "PARTIAL" in step_status or "WARNING" in step_status:
            badge    = '<span class="badge-warn">⚠️ PARTIAL</span>'
            row_bg   = "background:#fffbeb"
            extra    = ""
        else:
            badge    = '<span class="badge-pass">✅ PASS</span>'
            row_bg   = ""
            extra    = ""

        # ── Special biological result notes ──
        if key == "step6a_deseq2":
            de = s.get("significant_de_genes", 0)
            up = s.get("up_regulated", 0)
            dn = s.get("down_regulated", 0)
            extra += f'<br><small style="color:#475569">DE genes: {de:,} (↑{up:,} ↓{dn:,})</small>'
        elif key == "step7_pathway":
            kegg = s.get("enriched_kegg_pathways", 0)
            gobp = s.get("enriched_go_bp", 0)
            extra += f'<br><small style="color:#475569">KEGG: {kegg} pathways · GO-BP: {gobp} terms</small>'

        rows += f"""
        <tr style="{row_bg}">
            <td style="font-size:1.1rem;text-align:center">{icon}</td>
            <td style="text-align:center"><strong>{step_num}</strong></td>
            <td style="text-align:center">{step_name}</td>
            <td style="text-align:center"><code>{tool}</code></td>
            <td>{badge}{extra}</td>
            <td style="color:#475569">{tp_str}</td>
            <td style="color:#475569">{t:.3f}s</td>
        </tr>"""

    return rows



def _build_scrna_table(module_data):
    """Build scRNA-seq comparison table across Small/Medium/Large"""

    size_results = module_data.get("size_results", {})
    if not size_results:
        return "<tr><td colspan='9'>No data</td></tr>"

    step_labels = {
        "step1_cellranger" : ("🔬", "Cell Ranger",        "Alignment + Cell Calling", "32GB / 32GB / 64GB"),
        "step2_qc"         : ("📊", "QC Filtering",       "Filter low quality cells", "4GB / 16GB / 48GB"),
        "step3_norm"       : ("🔢", "Normalization + HVG","LogNormalize + ScaleData", "2GB / 16GB / 64GB"),
        "step4_dimred"     : ("📉", "PCA + UMAP",         "Dimensionality Reduction", "2GB / 8GB / 32GB"),
        "step5_clustering" : ("🔵", "Leiden Clustering",  "Community Detection",      "Low / Low / Medium"),
        "step6_markers"    : ("🏷️",  "Marker Genes",       "Wilcoxon Test per cluster","2GB / 8GB / 32GB"),
        "step7_trajectory" : ("🛤️",  "Trajectory",         "Monocle3 pseudotime",     "Low / Low / Medium"),
    }

    size_keys   = ["small", "medium", "large"]
    size_labels = ["Small (5k cells)", "Medium (30k cells)", "Large (100k cells)"]
    size_descs  = ["~5,000 cells", "~30,000 cells", "~100,000 cells"]
    ram_needed  = ["32GB", "32GB", "64GB"]

    # Header with 3 size columns
    rows = f"""
    <tr style="background:#f1f5f9">
        <th style="padding:8px;color:#475569">Icon</th>
        <th style="padding:8px;color:#475569">Step</th>
        <th style="padding:8px;color:#475569">Process</th>
        <th style="padding:8px;color:#475569">RAM Required</th>
        <th style="padding:8px;color:#475569;font-weight:600;text-align:center">
            Small Study<br><small style="font-weight:400;color:#64748b">{size_descs[0]}</small>
        </th>
        <th style="padding:8px;color:#475569;font-weight:600;text-align:center">
            Medium Study<br><small style="font-weight:400;color:#64748b">{size_descs[1]}</small>
        </th>
        <th style="padding:8px;background:#f1f5f9;color:#475569">
            Large Study<br><small style="font-weight:400;color:#64748b">{size_descs[2]}</small>
        </th>
    </tr>"""

    for key, (icon, step_name, process, ram_req) in step_labels.items():
        row = f"""
        <tr>
            <td style="font-size:1.1rem;padding:8px;text-align:center">{icon}</td>
            <td style="padding:8px"><strong>{step_name}</strong></td>
            <td style="padding:8px;color:#475569;font-size:0.8rem">{process}</td>
            <td style="padding:8px;font-size:0.75rem;color:#64748b">{ram_req}</td>"""

        for size_key in size_keys:
            if size_key not in size_results:
                row += "<td style='padding:8px'>N/A</td>"
                continue

            steps = size_results[size_key].get("pipeline_steps", {})
            if key not in steps:
                row += "<td style='padding:8px'>N/A</td>"
                continue

            s      = steps[key]
            status = str(s.get("status",""))
            note   = s.get("note","")

            if "FAIL" in status:
                cell_bg  = "background:#fff5f5"
                badge    = '<span class="badge-fail">🔴 FAIL</span>'
                note_str = f'<br><small style="color:#dc2626;font-size:0.7rem">{note[:60]}</small>' if note else ""
            else:
                cell_bg  = "background:#f0fdf4"
                badge    = '<span class="badge-pass">✅ PASS</span>'
                note_str = ""

                # Add biological info for PASS steps
                if key == "step2_qc":
                    surv = s.get("survival_rate", 0)
                    note_str = f'<br><small style="color:#475569;font-size:0.7rem">{surv}% survival</small>'
                elif key == "step5_clustering":
                    nc = s.get("n_clusters", 0)
                    note_str = f'<br><small style="color:#475569;font-size:0.7rem">{nc} clusters</small>'
                elif key == "step6_markers":
                    nm = s.get("total_markers", 0)
                    note_str = f'<br><small style="color:#475569;font-size:0.7rem">{nm:,} markers</small>'
                elif key == "step7_trajectory":
                    tg = s.get("trajectory_genes", 0)
                    note_str = f'<br><small style="color:#475569;font-size:0.7rem">{tg:,} traj. genes</small>'
                elif key == "step4_dimred":
                    un = s.get("umap_note","")
                    note_str = f'<br><small style="color:#475569;font-size:0.7rem">{un[:50]}</small>' if un else ""

            row += f'<td style="padding:8px;{cell_bg}">{badge}{note_str}</td>'

        row += "</tr>"
        rows += row

    # Score summary row
    rows += """<tr style="background:#f8fafc;border-top:2px solid #e2e8f0">
        <td colspan="4" style="padding:8px;font-weight:700;color:#1e293b">
            ⭐ Overall Score
        </td>"""

    for size_key in size_keys:
        if size_key in size_results:
            score = size_results[size_key].get("score", "N/A")
            cap   = size_results[size_key].get("capability","")
            color = "#dc2626" if score < 30 else "#d97706" if score < 60 else "#16a34a"
            rows += f"""<td style="padding:8px;font-weight:700;color:{color};font-size:1.1rem">
                {score}/100<br>
                <small style="font-size:0.7rem;font-weight:400;color:#64748b">{cap[:40]}</small>
            </td>"""
        else:
            rows += "<td>N/A</td>"

    rows += "</tr>"

    # RAM requirement summary
    rows += f"""<tr style="background:#eff6ff">
        <td colspan="4" style="padding:8px;font-weight:600;color:#1e40af">
            🧠 Your RAM: {module_data.get("hardware",{}).get("total_ram_gb","?")}GB
        </td>
        <td style="padding:8px;background:#f1f5f9;color:#475569;font-size:0.8rem">
            Needs 32GB<br>🔴 Insufficient
        </td>
        <td style="padding:8px;background:#f1f5f9;color:#475569;font-size:0.8rem">
            Needs 32GB<br>🔴 Insufficient
        </td>
        <td style="padding:8px;background:#f1f5f9;color:#475569;font-size:0.8rem">
            Needs 64GB<br>🔴 Insufficient
        </td>
    </tr>"""

    return rows



def _build_scrna_section(module_data, score_color):
    """Build complete scRNA-seq section with chart + comparison table"""

    size_results = module_data.get("size_results", {})
    hardware     = module_data.get("hardware", {})
    ram_gb       = hardware.get("total_ram_gb", "?")

    # Scores for each size
    small_score  = size_results.get("small",  {}).get("score", 0)
    medium_score = size_results.get("medium", {}).get("score", 0)
    large_score  = size_results.get("large",  {}).get("score", 0)
    overall      = module_data.get("score", 0)

    small_cap  = size_results.get("small",  {}).get("capability","")
    medium_cap = size_results.get("medium", {}).get("capability","")
    large_cap  = size_results.get("large",  {}).get("capability","")

    def sc(s):
        if not isinstance(s,(int,float)): return "#6b7280"
        return "#16a34a" if s>=80 else "#d97706" if s>=50 else "#dc2626"

    # Build comparison table rows
    step_labels = {
        "step1_cellranger" : ("🔬", "Cell Ranger",      "Alignment + Cell Calling", "32GB / 32GB / 64GB"),
        "step2_qc"         : ("📊", "QC Filtering",     "Filter low quality cells", "4GB / 16GB / 48GB"),
        "step3_norm"       : ("🔢", "Normalization+HVG","LogNormalize + ScaleData", "2GB / 16GB / 64GB"),
        "step4_dimred"     : ("📉", "PCA + UMAP",       "Dimensionality Reduction", "2GB / 8GB / 32GB"),
        "step5_clustering" : ("🔵", "Clustering",       "Leiden algorithm",         "Low / Low / Medium"),
        "step6_markers"    : ("🏷️",  "Marker Genes",     "Wilcoxon per cluster",     "2GB / 8GB / 32GB"),
        "step7_trajectory" : ("🛤️",  "Trajectory",       "Monocle3 pseudotime",      "Low / Low / Medium"),
    }

    size_keys = ["small", "medium", "large"]

    table_rows = ""
    for key, (icon, step_name, process, ram_req) in step_labels.items():
        row = f"""<tr>
            <td style="padding:8px;text-align:center;font-size:1.1rem">{icon}</td>
            <td style="padding:8px"><strong>{step_name}</strong></td>
            <td style="padding:8px;color:#475569;font-size:0.8rem">{process}</td>
            <td style="padding:8px;font-size:0.75rem;color:#64748b">{ram_req}</td>"""

        for sk in size_keys:
            steps = size_results.get(sk, {}).get("pipeline_steps", {})
            s     = steps.get(key, {})
            status = str(s.get("status",""))
            note   = s.get("note","")

            if "FAIL" in status:
                bg    = "background:#fff5f5"  # RED background for FAIL
                badge = '<span class="badge-fail" style="background:#fef2f2;color:#dc2626;padding:2px 8px;border-radius:99px;font-size:0.72rem;font-weight:600">🔴 FAIL</span>'
                short = note[:55] + "..." if len(note)>55 else note
                extra = f'<br><small style="color:#dc2626;font-size:0.7rem">{short}</small>' if note else ""
            else:
                bg    = "background:#f0fdf4"  # GREEN background for PASS
                badge = '<span class="badge-pass" style="background:#dcfce7;color:#166534;padding:2px 8px;border-radius:99px;font-size:0.72rem;font-weight:600">✅ PASS</span>'
                extra = ""
                if key == "step2_qc":
                    surv  = s.get("survival_rate",0)
                    cells = s.get("cells_passed",0)
                    extra = f'<br><small style="color:#475569;font-size:0.7rem">{cells:,} cells ({surv}%)</small>'
                elif key == "step4_dimred":
                    cells = s.get("cells_processed",0)
                    if cells > 0:
                        extra = f'<br><small style="color:#475569;font-size:0.7rem">{cells:,} cells processed</small>'
                elif key == "step5_clustering":
                    nc = s.get("n_clusters",0)
                    if nc > 0:
                        extra = f'<br><small style="color:#475569;font-size:0.7rem">{nc} clusters found</small>'
                elif key == "step6_markers":
                    nm = s.get("total_markers",0)
                    if nm > 0:
                        extra = f'<br><small style="color:#475569;font-size:0.7rem">{nm:,} markers</small>'
                elif key == "step7_trajectory":
                    tg = s.get("trajectory_genes",0)
                    if tg > 0:
                        extra = f'<br><small style="color:#475569;font-size:0.7rem">{tg:,} traj. genes</small>'

            row += f'<td style="padding:8px;{bg};text-align:center">{badge}{extra}</td>'

        row += "</tr>"
        table_rows += row

    # Score summary row
    score_row = """<tr style="background:#f8fafc;border-top:2px solid #e2e8f0">
        <td colspan="4" style="padding:10px;font-weight:700;font-size:0.9rem;text-align:center">
            ⭐ Overall Score
        </td>"""
    for sk, s_score, s_cap in [
        ("small", small_score, small_cap),
        ("medium", medium_score, medium_cap),
        ("large", large_score, large_cap)
    ]:
        color = sc(s_score)
        short_cap = s_cap[:45] + "..." if len(s_cap)>45 else s_cap
        score_row += f"""<td style="padding:10px;font-weight:800;
                            font-size:1.2rem;color:{color};text-align:center">
            {s_score}/100
            <br><small style="font-size:0.7rem;font-weight:400;
                              color:#64748b">{short_cap}</small>
        </td>"""
    score_row += "</tr>"

    # RAM row
    ram_row = f"""<tr style="background:#eff6ff">
        <td colspan="4" style="padding:8px;font-weight:600;color:#1e40af;text-align:center">
            🧠 Your RAM: {ram_gb}GB
        </td>
        <td style="padding:8px;background:#f1f5f9;color:#475569;font-size:0.8rem">
            Needs 32GB<br>🔴 Insufficient
        </td>
        <td style="padding:8px;background:#f1f5f9;color:#475569;font-size:0.8rem">
            Needs 32GB<br>🔴 Insufficient
        </td>
        <td style="padding:8px;background:#f1f5f9;color:#475569;font-size:0.8rem">
            Needs 64GB<br>🔴 Insufficient
        </td>
    </tr>"""

    # Unique chart ID
    chart_id = "scrnaChart"

    return f"""
    <div class="module-section">

        <!-- Header -->
        <div class="module-header">
            <div>
                <span style="font-size:1.5rem">🔬</span>
                <strong style="font-size:1.1rem;margin-left:0.5rem">
                    scRNA-seq Pipeline — Small / Medium / Large Dataset Comparison
                </strong>
            </div>
            <div style="text-align:right">
                <div style="font-size:1.8rem;font-weight:800;color:{score_color}">
                    {overall}/100
                </div>
                <div style="font-size:0.8rem;color:#64748b">Overall score</div>
            </div>
        </div>

        <!-- Score Chart -->
        <div style="display:grid;grid-template-columns:1fr 1fr;gap:1rem;margin-bottom:1.5rem">
            <div style="background:#f8fafc;border-radius:12px;padding:1.2rem;
                        border:1px solid #e2e8f0">
                <div style="font-size:0.85rem;font-weight:600;color:#475569;
                             margin-bottom:1rem">
                    📊 Score by Dataset Size
                </div>
                <canvas id="{chart_id}" height="120"></canvas>
            </div>
            <div style="background:#f8fafc;border-radius:12px;padding:1.2rem;
                        border:1px solid #e2e8f0">
                <div style="font-size:0.85rem;font-weight:600;color:#475569;
                             margin-bottom:0.75rem">
                    🧠 RAM Requirements vs Your Machine
                </div>
                <div style="display:flex;flex-direction:column;gap:0.6rem">
                    <div style="display:flex;justify-content:space-between;
                                align-items:center;font-size:0.85rem">
                        <span>Small Study (5k cells)</span>
                        <span style="font-weight:700;color:#dc2626">
                            Needs 32GB · You have {ram_gb}GB
                        </span>
                    </div>
                    <div style="height:6px;background:#e2e8f0;border-radius:3px">
                        <div style="width:{min(100,float(ram_gb)/32*100):.0f}%;
                                    height:100%;background:#dc2626;border-radius:3px">
                        </div>
                    </div>
                    <div style="display:flex;justify-content:space-between;
                                align-items:center;font-size:0.85rem">
                        <span>Medium Study (30k cells)</span>
                        <span style="font-weight:700;color:#dc2626">
                            Needs 32GB · You have {ram_gb}GB
                        </span>
                    </div>
                    <div style="height:6px;background:#e2e8f0;border-radius:3px">
                        <div style="width:{min(100,float(ram_gb)/32*100):.0f}%;
                                    height:100%;background:#dc2626;border-radius:3px">
                        </div>
                    </div>
                    <div style="display:flex;justify-content:space-between;
                                align-items:center;font-size:0.85rem">
                        <span>Large Study (100k cells)</span>
                        <span style="font-weight:700;color:#dc2626">
                            Needs 64GB · You have {ram_gb}GB
                        </span>
                    </div>
                    <div style="height:6px;background:#e2e8f0;border-radius:3px">
                        <div style="width:{min(100,float(ram_gb)/64*100):.0f}%;
                                    height:100%;background:#dc2626;border-radius:3px">
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Comparison Table -->
        <div class="table-wrap">
            <table>
                <thead>
                    <tr style="background:#f1f5f9">
                        <th style="padding:8px">Icon</th>
                        <th style="padding:8px">Step</th>
                        <th style="padding:8px">Process</th>
                        <th style="padding:8px">RAM Required</th>
                        <th style="padding:8px;background:#f1f5f9;color:#475569;font-weight:600">
                            Small Study<br>
                            <span style="font-weight:400;color:#64748b;font-size:0.75rem">~5,000 cells</span>
                        </th>
                        <th style="padding:8px;background:#f1f5f9;color:#475569;font-weight:600">
                            Medium Study<br>
                            <span style="font-weight:400;color:#64748b;font-size:0.75rem">~30,000 cells</span>
                        </th>
                        <th style="padding:8px;background:#f1f5f9;color:#475569;font-weight:600">
                            Large Study<br>
                            <span style="font-weight:400;color:#64748b;font-size:0.75rem">~100,000 cells</span>
                        </th>
                    </tr>
                </thead>
                <tbody>
                    {table_rows}
                    {score_row}
                    {ram_row}
                </tbody>
            </table>
        </div>

    </div>

    <script>
    new Chart(document.getElementById('{chart_id}'), {{
        type: 'bar',
        data: {{
            labels: ['Small (5k)', 'Medium (30k)', 'Large (100k)'],
            datasets: [{{
                label: 'Score',
                data: [{small_score}, {medium_score}, {large_score}],
                backgroundColor: [
                    '{sc(small_score)}',
                    '{sc(medium_score)}',
                    '{sc(large_score)}'
                ],
                borderRadius: 8,
            }}]
        }},
        options: {{
            responsive: true,
            scales: {{
                y: {{ beginAtZero: true, max: 100 }}
            }},
            plugins: {{ legend: {{ display: false }} }}
        }}
    }});
    </script>
    """





def _build_generic_table(module_data, step_labels):
    """Generic pipeline table builder for any module"""
    rows  = ""
    steps = module_data.get("pipeline_steps", {})

    for key, (icon, step_num, step_name, tool) in step_labels.items():
        if key not in steps:
            continue
        s      = steps[key]
        t      = s.get("time_seconds", 0)
        tp     = s.get("throughput_reads_per_sec", 0)
        tp_str = f"{tp:,}/s" if tp else "N/A"
        status = str(s.get("status",""))
        note   = s.get("note","")

        if "FAIL" in status:
            badge  = '<span class="badge-fail">🔴 FAIL</span>'
            row_bg = "background:#fff5f5"
            extra  = f'<br><small style="color:#dc2626">{note[:70]}</small>' if note else ""
        else:
            badge  = '<span class="badge-pass">✅ PASS</span>'
            row_bg = ""
            extra  = f'<br><small style="color:#475569">{note[:70]}</small>' if note else ""

        rows += f"""
        <tr style="{row_bg}">
            <td style="font-size:1.1rem;text-align:center">{icon}</td>
            <td><strong>{step_num}</strong></td>
            <td>{step_name}</td>
            <td><code>{tool}</code></td>
            <td>{badge}{extra}</td>
            <td style="color:#475569">{tp_str}</td>
            <td style="color:#475569">{t:.3f}s</td>
        </tr>"""

    return rows


def _build_epigenomics_table(module_data):
    """Build Epigenomics pipeline table — ChIP-seq + ATAC-seq"""
    rows  = ""
    steps = module_data.get("pipeline_steps", {})

    step_labels = {
        "chipseq_qc"      : ("📊", "ChIP-1", "QC",                   "FastQC"),
        "chipseq_align"   : ("🔗", "ChIP-2", "Alignment",            "Bowtie2"),
        "chipseq_peaks"   : ("🏔️",  "ChIP-3", "Peak Calling",         "MACS3"),
        "chipseq_motif"   : ("🔤", "ChIP-4", "Motif Analysis",       "HOMER"),
        "chipseq_diff"    : ("📈", "ChIP-5", "Differential Binding", "DiffBind"),
        "atacseq_qc"      : ("📊", "ATAC-1", "QC + Fragment Size",   "FastQC"),
        "atacseq_align"   : ("🔗", "ATAC-2", "Alignment",            "Bowtie2"),
        "atacseq_peaks"   : ("🏔️",  "ATAC-3", "Peak Calling",         "MACS3"),
        "atacseq_chromvar": ("🔍", "ATAC-4", "TF Activity",          "chromVAR"),
        "deeptools"       : ("📊", "BOTH-5", "Coverage + Viz",       "deepTools"),
    }

    for key, (icon, step_num, step_name, tool) in step_labels.items():
        if key not in steps:
            continue
        s      = steps[key]
        t      = s.get("time_seconds", 0)
        tp     = s.get("throughput_reads_per_sec", 0)
        tp_str = f"{tp:,}/s" if tp else "N/A"
        status = str(s.get("status",""))
        note   = s.get("note","")

        if "FAIL" in status:
            badge  = '<span class="badge-fail">🔴 FAIL</span>'
            row_bg = "background:#fff5f5"
            extra  = f'<br><small style="color:#dc2626">{note[:60]}</small>' if note else ""
        else:
            badge  = '<span class="badge-pass">✅ PASS</span>'
            row_bg = ""
            extra  = ""

            # Add biological highlights
            if key == "chipseq_peaks":
                n = s.get("peaks_called", 0)
                extra = f'<br><small style="color:#475569">{n:,} peaks called</small>'
            elif key == "chipseq_motif":
                n = s.get("motifs_found", 0)
                extra = f'<br><small style="color:#475569">{n} motifs found</small>'
            elif key == "chipseq_diff":
                n = s.get("differential_peaks", 0)
                extra = f'<br><small style="color:#475569">{n:,} differential peaks</small>'
            elif key == "atacseq_peaks":
                n   = s.get("peaks_called", 0)
                tss = s.get("tss_enrichment", 0)
                extra = f'<br><small style="color:#475569">{n:,} peaks · TSS enrichment: {tss}</small>'
            elif key == "atacseq_chromvar":
                n = s.get("significant_tfs", 0)
                extra = f'<br><small style="color:#475569">{n} significant TF activities</small>'
            elif key == "atacseq_qc":
                sub = s.get("sub_nucleosomal_pct", 0)
                extra = f'<br><small style="color:#475569">Sub-nucleosomal: {sub}% · Nucleosomal pattern: ✅</small>'

        rows += f"""
        <tr style="{row_bg}">
            <td style="font-size:1.1rem;text-align:center">{icon}</td>
            <td><strong>{step_num}</strong></td>
            <td>{step_name}</td>
            <td><code>{tool}</code></td>
            <td>{badge}{extra}</td>
            <td style="color:#475569">{tp_str}</td>
            <td style="color:#475569">{t:.3f}s</td>
        </tr>"""

    return rows


def _build_protein_table(module_data):
    """Build detailed Protein pipeline table"""
    rows  = ""
    steps = module_data.get("pipeline_steps", {})

    step_labels = {
        "step1_seq_prep" : ("🧬", "Step 1", "Sequence Preparation", "FASTA parsing"),
        "step2_msa"      : ("🔍", "Step 2", "MSA Search",           "HHblits / JackHMMER"),
        "step3_alphafold": ("🔮", "Step 3", "Structure Prediction",  "AlphaFold2"),
        "step4_esmfold"  : ("⚡", "Step 4", "Structure Prediction",  "ESMFold (Meta AI)"),
        "step5_variant"  : ("🧪", "Step 5", "Variant Effect",        "FoldX / DynaMut2"),
        "step6_blast"    : ("🔎", "Step 6", "Sequence Search",       "BLAST"),
    }

    for key, (icon, step_num, step_name, tool) in step_labels.items():
        if key not in steps:
            continue
        s          = steps[key]
        t          = s.get("time_seconds", 0)
        tp         = s.get("throughput_reads_per_sec", 0)
        tp_str     = f"{tp:,}/s" if tp else "N/A"
        status     = str(s.get("status",""))
        note       = s.get("note","")

        if "FAIL" in status:
            badge  = '<span class="badge-fail">🔴 FAIL</span>'
            row_bg = "background:#fff5f5"
            extra  = f'<br><small style="color:#dc2626">{note[:60]}</small>' if note else ""
        else:
            badge  = '<span class="badge-pass">✅ PASS</span>'
            row_bg = ""
            extra  = f'<br><small style="color:#475569">{note[:60]}</small>' if note else ""

        # Add specific biological info
        if key == "step3_alphafold" and "FAIL" not in status:
            plddt = s.get("mean_plddt_score", 0)
            high  = s.get("high_confidence", 0)
            extra += f'<br><small style="color:#475569">Mean pLDDT: {plddt} · High conf: {high} proteins</small>'
        elif key == "step4_esmfold" and "FAIL" not in status:
            plddt = s.get("mean_plddt_score", 0)
            extra += f'<br><small style="color:#475569">Mean pLDDT: {plddt} · No MSA needed</small>'
        elif key == "step5_variant" and "FAIL" not in status:
            path  = s.get("predicted_pathogenic", 0)
            total = s.get("variants_analyzed", 0)
            extra += f'<br><small style="color:#475569">Predicted pathogenic: {path}/{total} variants</small>'
        elif key == "step2_msa" and "FAIL" not in status:
            depth = s.get("mean_msa_depth", 0)
            extra += f'<br><small style="color:#475569">Mean MSA depth: {depth:,} sequences</small>'

        rows += f"""
        <tr style="{row_bg}">
            <td style="font-size:1.1rem;text-align:center">{icon}</td>
            <td><strong>{step_num}</strong></td>
            <td>{step_name}</td>
            <td><code>{tool}</code></td>
            <td>{badge}{extra}</td>
            <td style="color:#475569">{tp_str}</td>
            <td style="color:#475569">{t:.3f}s</td>
        </tr>"""

    return rows


def generate_html_report(system_profile, benchmark_results,
                         output_dir="results"):
    os.makedirs(output_dir, exist_ok=True)
    timestamp   = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(
        output_dir, f"biomark_report_{timestamp}.html"
    )

    # ── System Info ──
    cpu     = system_profile.get("cpu", {})
    ram     = system_profile.get("ram", {})
    storage = system_profile.get("storage", {})
    gpu     = system_profile.get("gpu", {})
    os_info = system_profile.get("os", {})

    chip      = cpu.get("apple_chip", cpu.get("brand","Unknown"))
    ram_total = ram.get("total_gb","N/A")
    ram_used  = ram.get("used_gb","N/A")
    ram_pct   = ram.get("percent_used", 0)
    ssd_total = storage.get("total_gb","N/A")
    ssd_free  = storage.get("free_gb","N/A")
    ssd_pct   = storage.get("percent_used", 0)
    gpu_model = gpu.get("model","N/A")
    macos     = os_info.get("macos_version","N/A")
    python_v  = os_info.get("python_version","N/A")
    mac_model = cpu.get("mac_model","N/A")
    cores     = cpu.get("physical_cores","N/A")

    # ── Dynamic Time Estimates ──
    time_estimates = {
        "DNA_WES"     : [
            ("🧬", "WES single sample (BWA-MEM2 + GATK)", "4–8 hours",   "#185FA5", 80),
            ("🧬", "WES trio analysis (3 samples)",        "12–24 hours", "#d97706", 95),
        ],
        "RNA_seq"     : [
            ("🔬", "Bulk RNA-seq 6 samples (STAR + DESeq2)", "3–6 hours", "#185FA5", 60),
        ],
        "scRNA_seq"   : [
            ("🔬", "scRNA-seq small 5k cells (Cell Ranger)", "2–4 hours",  "#185FA5", 40),
            ("🔬", "scRNA-seq medium 30k cells",             "6–12 hours", "#d97706", 70),
            ("🔬", "scRNA-seq large 100k cells",             "24–48 hours","#dc2626", 100),
        ],
        "Protein"     : [
            ("🔮", "AlphaFold2 per protein (CPU only)", "30 min–24 hrs", "#d97706", 70),
            ("⚡", "ESMFold per protein",               "1–10 min",      "#185FA5", 20),
        ],
        "Epigenomics" : [
            ("🧬", "ChIP-seq single sample (Bowtie2+MACS3)", "2–4 hours",  "#185FA5", 40),
            ("🧬", "ATAC-seq single sample",                 "3–5 hours",  "#185FA5", 50),
            ("🧬", "Differential binding (DiffBind)",        "1–3 hours",  "#d97706", 30),
        ],
        "Metagenomics": [
            ("🦠", "Kraken2 classification (per sample)",    "30min–2hrs", "#185FA5", 30),
            ("🌿", "QIIME2 16S analysis",                   "1–4 hours",  "#185FA5", 40),
            ("🔬", "HUMAnN3 functional profiling",          "2–8 hours",  "#d97706", 60),
        ],
        "LongRead"    : [
            ("📡", "Nanopore basecalling (Dorado GPU)",      "1–6 hours",  "#185FA5", 50),
            ("🔗", "Minimap2 alignment",                     "30min–4hrs", "#185FA5", 35),
            ("🧬", "Sniffles2 SV calling",                   "1–3 hours",  "#d97706", 30),
        ],
        "Assembly"    : [
            ("🔧", "SPAdes short read assembly",             "2–48 hours", "#d97706", 70),
            ("🔧", "Flye long read assembly",                "4–72 hours", "#dc2626", 90),
            ("🔧", "Hifiasm HiFi assembly",                  "6–48 hours", "#dc2626", 80),
        ],
    }

    time_estimates_html = ""
    for module_name in benchmark_results.keys():
        for key, estimates in time_estimates.items():
            # Exact match only to avoid duplicates
            if key.lower() == module_name.lower():
                for icon, label, duration, color, width in estimates:
                    time_estimates_html += f"""
                <div style="display:flex;justify-content:space-between;align-items:center">
                    <span>{icon} {label}</span>
                    <span style="font-weight:600;color:{color}">{duration}</span>
                </div>
                <div style="height:4px;background:#e2e8f0;border-radius:2px">
                    <div style="width:{width}%;height:100%;background:{color};border-radius:2px"></div>
                </div>"""

    if not time_estimates_html:
        time_estimates_html = '<div style="color:#94a3b8;font-size:0.85rem">Run benchmarks to see time estimates</div>'

    # ── Overall Score ──
    scores = [
        d.get("score",0)
        for d in benchmark_results.values()
        if isinstance(d.get("score"), (int,float))
    ]
    overall = round(sum(scores)/len(scores), 1) if scores else 0

    # ── Score Color ──
    overall_color = _score_color(overall)

    # ── Hardware Warnings Section ──
    hw_warnings_html = ""
    for module_data in benchmark_results.values():
        for w in module_data.get("hardware_warnings", []):
            level  = w.get("level","")
            comp   = w.get("component","")
            msg    = w.get("message","")
            rec    = w.get("recommendation","")
            color  = (
                "#fef2f2" if "🔴" in level else
                "#fffbeb" if "🟡" in level else
                "#f0fdf4"
            )
            border = (
                "#fca5a5" if "🔴" in level else
                "#fcd34d" if "🟡" in level else
                "#86efac"
            )
            hw_warnings_html += f"""
            <div style="background:{color};border-left:4px solid {border};
                        padding:0.75rem 1rem;border-radius:0 8px 8px 0;
                        margin-bottom:0.5rem;">
                <strong>{level} [{comp}]</strong> {msg}<br>
                <small style="color:#6b7280">💡 {rec}</small>
            </div>"""
        break  # Only show once

    # ── Module Sections ──
    module_sections_html = ""
    for module_name, module_data in benchmark_results.items():
        score      = module_data.get("score", "N/A")
        capability = module_data.get("capability","")
        total_time = module_data.get("total_time_seconds","N/A")
        sc         = _score_color(score)
        sb         = _score_badge(score)

        # Build table based on module type
        if "DNA" in module_name or "WES" in module_name:
            table_headers = """
                <th>Icon</th><th>Step</th><th>Process</th>
                <th>Tool</th><th>Status</th><th>Bottleneck</th>
                <th>Throughput</th><th>Time</th>"""
            table_rows = _build_wes_table(module_data)
            module_icon = "🧬"
            module_label = "WES Pipeline — Rare Disease Analysis"
        elif "scRNA" in module_name:
            table_headers = ""
            table_rows    = _build_scrna_table(module_data)
            module_icon   = "🔬"
            module_label  = "scRNA-seq Pipeline — Small / Medium / Large Study Comparison"
        elif "RNA" in module_name:
            table_headers = """
                <th>Icon</th><th>Step</th><th>Process</th>
                <th>Tool</th><th>Status / Notes</th>
                <th>Throughput</th><th>Time</th>"""
            table_rows  = _build_rna_table(module_data)
            module_icon = "🔬"
            module_label = "Bulk RNA-seq Pipeline — 3 Disease vs 3 Control"
        elif "Epigenomics" in module_name:
            table_headers = """
                <th>Icon</th><th>Step</th><th>Process</th>
                <th>Tool</th><th>Status / Notes</th>
                <th>Throughput</th><th>Time</th>"""
            table_rows  = _build_epigenomics_table(module_data)
            module_icon = "🧬"
            module_label = "Epigenomics — ChIP-seq + ATAC-seq Pipeline"
        elif "Protein" in module_name:
            table_headers = """
                <th>Icon</th><th>Step</th><th>Process</th>
                <th>Tool</th><th>Status / Notes</th>
                <th>Throughput</th><th>Time</th>"""
            table_rows  = _build_protein_table(module_data)
            module_icon = "🔮"
            module_label = "Protein Structure — AlphaFold2 + ESMFold + BLAST"
        elif "Metagenomics" in module_name:
            table_headers = """
                <th>Icon</th><th>Step</th><th>Process</th>
                <th>Tool</th><th>Status / Notes</th>
                <th>Throughput</th><th>Time</th>"""
            table_rows  = _build_generic_table(module_data, {
                "step1_qc"       : ("📊","Step 1","QC","FastQC + fastp"),
                "step2_kraken2"  : ("🦠","Step 2","Taxonomic Classification","Kraken2"),
                "step3_bracken"  : ("📊","Step 3","Abundance Estimation","Bracken"),
                "step4_metaphlan": ("🧬","Step 4","Microbial Profiling","MetaPhlAn4"),
                "step5_qiime2"   : ("🌿","Step 5","16S rRNA Analysis","QIIME2"),
                "step6_humann"   : ("🔬","Step 6","Functional Profiling","HUMAnN3"),
            })
            module_icon  = "🦠"
            module_label = "Metagenomics — Kraken2 + MetaPhlAn4 + QIIME2 + HUMAnN3"

        elif "LongRead" in module_name:
            table_headers = """
                <th>Icon</th><th>Step</th><th>Process</th>
                <th>Tool</th><th>Status / Notes</th>
                <th>Throughput</th><th>Time</th>"""
            table_rows  = _build_generic_table(module_data, {
                "step1_basecalling": ("📡","Step 1","Basecalling","Dorado (Nanopore)"),
                "step2_qc"         : ("📊","Step 2","QC","NanoPlot + FastQC"),
                "step3_minimap2"   : ("🔗","Step 3","Alignment","Minimap2"),
                "step4_sv_calling" : ("🧬","Step 4","SV Calling","Sniffles2"),
                "step5_medaka"     : ("🔧","Step 5","Variant + Polishing","Medaka"),
                "step6_methylation": ("🧪","Step 6","Methylation Detection","Modbam2bed"),
            })
            module_icon  = "🧬"
            module_label = "Long Read Sequencing — Nanopore + PacBio Pipeline"

        elif "Assembly" in module_name:
            table_headers = """
                <th>Icon</th><th>Step</th><th>Process</th>
                <th>Tool</th><th>Status / Notes</th>
                <th>Throughput</th><th>Time</th>"""
            table_rows  = _build_generic_table(module_data, {
                "step1_spades" : ("🔧","Step 1","Short Read Assembly","SPAdes"),
                "step2_flye"   : ("🔧","Step 2","Long Read Assembly","Flye"),
                "step3_hifiasm": ("🔧","Step 3","HiFi Assembly","Hifiasm"),
                "step4_quast"  : ("📊","Step 4","Assembly QC","QUAST"),
                "step5_busco"  : ("🔍","Step 5","Completeness Check","BUSCO"),
            })
            module_icon  = "🔨"
            module_label = "Genome Assembly — SPAdes + Flye + Hifiasm + BUSCO"

        else:
            table_headers = "<th>Module</th><th>Score</th>"
            table_rows    = ""
            module_icon   = "🔬"
            module_label  = module_name

        module_sections_html += f"""
        <div class="module-section">
            <div class="module-header">
                <div>
                    <span style="font-size:1.5rem">{module_icon}</span>
                    <strong style="font-size:1.1rem;margin-left:0.5rem">
                        {module_label}
                    </strong>
                </div>
                <div style="text-align:right">
                    <div style="font-size:1.8rem;font-weight:800;color:{sc}">
                        {score}/100
                    </div>
                    <div style="font-size:0.8rem;color:#64748b">
                        Total time: {total_time}s
                    </div>
                </div>
            </div>
            {"<div class='capability'>"+capability+"</div>" if capability else ""}
            <div class="table-wrap">
                <table>
                    <thead><tr>{table_headers}</tr></thead>
                    <tbody>{table_rows}</tbody>
                </table>
            </div>
        </div>"""

    # ── Chart Data — expand scRNA into 3 bars ──
    chart_labels = []
    chart_scores = []
    chart_times  = []
    chart_colors = []

    for module_name, d in benchmark_results.items():
        if "scRNA" in module_name:
            # Show only overall score in main chart
            # Size breakdown is already shown in the detailed table
            chart_labels.append("scRNA-seq")
            chart_scores.append(d.get("score", 0))
            chart_times.append(d.get("total_time_seconds", 0))
        else:
            chart_labels.append(module_name)
            chart_scores.append(d.get("score", 0))
            chart_times.append(d.get("total_time_seconds", 0))

    chart_labels = json.dumps(chart_labels)
    chart_scores = json.dumps(chart_scores)
    chart_times  = json.dumps(chart_times)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width,initial-scale=1.0">
    <title>BioMark Report — {datetime.now().strftime('%Y-%m-%d')}</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        *{{box-sizing:border-box;margin:0;padding:0}}
        body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;
              background:#f8fafc;color:#1e293b;line-height:1.6}}

        /* Header */
        .header{{background:linear-gradient(135deg,#1a3a5c,#185FA5);
                 color:white;padding:2.5rem 2rem;text-align:center}}
        .header h1{{font-size:2rem;font-weight:700}}
        .overall-score{{font-size:3.5rem;font-weight:800;
                        color:{overall_color};margin:0.5rem 0;
                        text-shadow:0 2px 4px rgba(0,0,0,0.3)}}
        .header .meta{{font-size:0.85rem;opacity:0.75;margin-top:0.5rem}}

        /* Container */
        .container{{max-width:1200px;margin:0 auto;padding:2rem 1.5rem}}

        /* Section title */
        .section-title{{font-size:1rem;font-weight:600;color:#1e293b;
                        margin:2rem 0 1rem;padding-bottom:0.4rem;
                        border-bottom:2px solid #e2e8f0;display:flex;
                        align-items:center;gap:0.5rem}}

        /* System cards */
        .cards{{display:grid;
                grid-template-columns:repeat(auto-fit,minmax(160px,1fr));
                gap:1rem;margin-bottom:1.5rem}}
        .card{{background:white;border-radius:12px;padding:1.1rem;
               box-shadow:0 1px 4px rgba(0,0,0,0.07);
               border:1px solid #e2e8f0}}
        .card-icon{{font-size:1.4rem;margin-bottom:0.4rem}}
        .card-label{{font-size:0.7rem;color:#64748b;text-transform:uppercase;
                     letter-spacing:0.05em;font-weight:600}}
        .card-value{{font-size:1rem;font-weight:700;color:#1e293b;
                     margin-top:0.2rem}}
        .card-sub{{font-size:0.75rem;color:#94a3b8;margin-top:0.15rem}}
        .progress-bar{{height:5px;background:#e2e8f0;border-radius:3px;
                       overflow:hidden;margin-top:0.4rem}}
        .progress-fill{{height:100%;border-radius:3px}}
        .fill-red{{background:#ef4444}}
        .fill-amber{{background:#f59e0b}}
        .fill-green{{background:#10b981}}

        /* Charts */
        .charts{{display:grid;grid-template-columns:1fr 1fr;
                 gap:1.5rem;margin-bottom:2rem}}
        .chart-card{{background:white;border-radius:12px;padding:1.5rem;
                     box-shadow:0 1px 4px rgba(0,0,0,0.07);
                     border:1px solid #e2e8f0}}
        .chart-title{{font-size:0.85rem;font-weight:600;color:#475569;
                      margin-bottom:1rem}}

        /* Module sections */
        .module-section{{background:white;border-radius:16px;
                          padding:1.5rem;margin-bottom:2rem;
                          box-shadow:0 2px 8px rgba(0,0,0,0.08);
                          border:1px solid #e2e8f0}}
        .module-header{{display:flex;justify-content:space-between;
                        align-items:center;margin-bottom:1rem;
                        padding-bottom:1rem;
                        border-bottom:1px solid #e2e8f0}}
        .capability{{background:#f8fafc;border-radius:8px;
                     padding:0.6rem 1rem;font-size:0.9rem;
                     font-weight:500;margin-bottom:1rem;
                     border-left:4px solid #185FA5}}

        /* Table */
        .table-wrap{{overflow-x:auto;border-radius:8px;
                     border:1px solid #e2e8f0}}
        table{{width:100%;border-collapse:collapse;font-size:0.8rem}}
        thead{{background:#f1f5f9}}
        th{{padding:0.6rem 0.75rem;text-align:left;font-weight:600;
            color:#475569;font-size:0.72rem;text-transform:uppercase;
            letter-spacing:0.04em;white-space:nowrap}}
        td{{padding:0.6rem 0.75rem;border-top:1px solid #f1f5f9;
            color:#334155;vertical-align:top}}
        tr:hover td{{background:#f8fafc}}
        code{{background:#f1f5f9;padding:2px 6px;border-radius:4px;
              font-size:0.75rem;color:#0369a1}}

        /* Badges */
        .badge-pass{{background:#dcfce7;color:#166534;padding:2px 8px;
                     border-radius:99px;font-size:0.72rem;font-weight:600;
                     white-space:nowrap}}
        .badge-fail{{background:#f1f5f9;color:#475569;padding:2px 8px;
                     border-radius:99px;font-size:0.72rem;font-weight:600;
                     white-space:nowrap}}
        .badge-warn{{background:#fef9c3;color:#854d0e;padding:2px 8px;
                     border-radius:99px;font-size:0.72rem;font-weight:600;
                     white-space:nowrap}}

        /* Footer */
        .footer{{text-align:center;padding:2rem;color:#94a3b8;
                 font-size:0.85rem}}
        .footer a{{color:#185FA5;text-decoration:none}}
    </style>
</head>
<body>

<!-- HEADER -->
<div class="header">
    <h1>🧬 BioMark Benchmark Report</h1>
    <p style="opacity:0.8">Bioinformatics Hardware Performance Analysis</p>
    <div class="overall-score">{overall}/100</div>
    <div style="font-size:0.9rem;opacity:0.85">Overall BioMark Score</div>
    <div class="meta">
        {datetime.now().strftime('%Y-%m-%d %H:%M')} &nbsp;|&nbsp;
        {mac_model} &nbsp;|&nbsp; macOS {macos}
    </div>
</div>

<div class="container">

    <!-- SYSTEM PROFILE -->
    <div class="section-title">🖥️ System Profile</div>
    <div class="cards">
        <div class="card">
            <div class="card-icon">⚡</div>
            <div class="card-label">Processor</div>
            <div class="card-value">{chip}</div>
            <div class="card-sub">{cores} CPU cores · arm64</div>
        </div>
        <div class="card">
            <div class="card-icon">🧠</div>
            <div class="card-label">Memory (RAM)</div>
            <div class="card-value">{ram_total} GB</div>
            <div class="card-sub">{ram_used}GB used · {ram_pct}%</div>
            <div class="progress-bar">
                <div class="progress-fill
                    {'fill-red' if ram_pct>80 else 'fill-amber' if ram_pct>60 else 'fill-green'}"
                    style="width:{ram_pct}%"></div>
            </div>
        </div>
        <div class="card">
            <div class="card-icon">💾</div>
            <div class="card-label">Storage (SSD)</div>
            <div class="card-value">{ssd_total} GB</div>
            <div class="card-sub">{ssd_free}GB free · {ssd_pct}% used</div>
            <div class="progress-bar">
                <div class="progress-fill
                    {'fill-red' if ssd_pct>80 else 'fill-amber' if ssd_pct>60 else 'fill-green'}"
                    style="width:{ssd_pct}%"></div>
        </div>
        </div>
        <div class="card">
            <div class="card-icon">🎮</div>
            <div class="card-label">GPU</div>
            <div class="card-value">{gpu_model}</div>
            <div class="card-sub">Unified Memory · Metal</div>
        </div>
        <div class="card">
            <div class="card-icon">🐍</div>
            <div class="card-label">Python</div>
            <div class="card-value">{python_v}</div>
            <div class="card-sub">macOS {macos}</div>
        </div>
    </div>

    <!-- HARDWARE WARNINGS -->
    <div class="section-title">⚠️ Hardware Assessment</div>
    {hw_warnings_html}

    <!-- CHARTS -->
    <div class="section-title">📊 Score Overview</div>
    <div class="charts">
        <div class="chart-card">
            <div class="chart-title">Module Scores (0–100) · Higher is better</div>
            <canvas id="scoreChart"></canvas>
        </div>
        <div class="chart-card">
            <div class="chart-title">⏱️ Real Pipeline Time Estimates (on recommended hardware)</div>
            <div style="font-size:0.8rem;color:#64748b;margin-bottom:0.75rem">
                Approximate wall-clock time on a machine with sufficient RAM
            </div>
            <div style="display:flex;flex-direction:column;gap:0.6rem;font-size:0.82rem">
                {time_estimates_html}
            </div>
        </div>
    </div>

    <!-- MODULE SECTIONS -->
    <div class="section-title">🔬 Detailed Pipeline Results</div>
    {module_sections_html}

</div>

<!-- FOOTER -->
<div class="footer">
    Generated by
    <a href="https://github.com/shahbazigenomics/BioMark">BioMark</a>
    — Open-source Bioinformatics Benchmark Tool by
    <strong>Amir Shahbazi</strong> · PhD Candidate in Genomics
</div>

<script>
const labels = {chart_labels};
const scores = {chart_scores};
const times  = {chart_times};

new Chart(document.getElementById('scoreChart'),{{
    type:'bar',
    data:{{
        labels:labels,
        datasets:[{{
            label:'Score',
            data:scores,
            backgroundColor:scores.map(s=>
                s===0?'#e2e8f0':s>=80?'#16a34a':s>=50?'#d97706':'#dc2626'
            ),
            borderRadius:8,
        }}]
    }},
    options:{{
        responsive:true,
        scales:{{y:{{beginAtZero:true,max:100}}}},
        plugins:{{
            legend:{{display:false}},
            tooltip:{{
                callbacks:{{
                    label: function(ctx) {{
                        return ctx.raw === 0
                            ? 'Cannot run on this machine'
                            : 'Score: ' + ctx.raw + '/100';
                    }}
                }}
            }}
        }}
    }}
}});

new Chart(document.getElementById('timeChart'),{{
    type:'bar',
    data:{{
        labels:labels,
        datasets:[{{
            label:'Time (s)',
            data:times.map(t => t > 0 ? t : null),
            backgroundColor:'#185FA5',
            borderRadius:8,
        }}]
    }},
    options:{{
        responsive:true,
        scales:{{
            y:{{
                beginAtZero:true,
                title:{{display:true, text:'seconds'}}
            }}
        }},
        plugins:{{
            legend:{{display:false}},
            tooltip:{{
                callbacks:{{
                    label: function(ctx) {{
                        return ctx.raw === null ? 'Not measured' : ctx.raw + 's';
                    }}
                }}
            }}
        }}
    }}
}});

// Add "Not measured" text labels on zero bars
Chart.register({{
    id: 'notMeasuredLabels',
    afterDraw(chart) {{
        const ctx2 = chart.ctx;
        chart.data.datasets.forEach((dataset, i) => {{
            const meta = chart.getDatasetMeta(i);
            meta.data.forEach((bar, index) => {{
                const originalVal = times[index];
                if (originalVal === 0) {{
                    ctx2.save();
                    ctx2.fillStyle = '#94a3b8';
                    ctx2.font = '11px sans-serif';
                    ctx2.textAlign = 'center';
                    ctx2.fillText('Not measured', bar.x, chart.chartArea.bottom - 10);
                    ctx2.restore();
                }}
            }});
        }});
    }}
}});
</script>
</body>
</html>"""

    with open(output_file, "w") as f:
        f.write(html)

    print(f"\n📄 HTML Report: {output_file}")
    return output_file
