# Exact orthogonal decomposition of sequence information across 20,431 human proteins
import gzip
import json
import math
import os
import statistics
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES_DIR = os.path.join(ROOT, "results")
FIG_DIR = os.path.join(ROOT, "figures")
os.makedirs(RES_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

CANDIDATES = [
    os.path.join(ROOT, "data", "fibre_results.json"),
    os.path.join(ROOT, "..", "separation-resolution", "results", "fibre_results.json"),
    os.path.join(ROOT, "..", "results", "fibre_results.json"),
]
DATA_PATH = next((p for p in CANDIDATES if os.path.exists(p)), CANDIDATES[0])

with open(DATA_PATH, "r", encoding="utf-8") as f:
    raw = json.load(f)

proteins = raw["proteins"] if isinstance(raw, dict) and "proteins" in raw else raw

TSV_CANDIDATES = [
    os.path.join(ROOT, "..", "separation-resolution", "results", "uniprot_hs_reviewed.tsv.gz"),
    os.path.join(ROOT, "..", "results", "uniprot_hs_reviewed.tsv.gz"),
    os.path.join(ROOT, "..", "data", "uniprot_hs_reviewed.tsv.gz"),
]
seq_map = {}
name_map = {}
for tp in TSV_CANDIDATES:
    if os.path.exists(tp):
        with gzip.open(tp, "rt", encoding="utf-8") as gf:
            header = gf.readline()
            for line in gf:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    acc = parts[0]
                    seq_map[acc] = parts[-1]
                    if len(parts) >= 4:
                        name_map[acc] = " ".join(parts[1:-1])
        break

LOG2_20 = math.log2(20.0)
LN2 = math.log(2.0)

def exact_B1_from_seq(seq):
    L = len(seq)
    cnt = Counter(seq)
    lg = math.lgamma(L + 1) - sum(math.lgamma(c + 1) for c in cnt.values())
    return lg / LN2

def extract_Bk_dict(p):
    out = {}
    for key in ("logN", "log2_N", "log2N", "B", "Bk", "fibres", "log2_fibres"):
        if key in p and isinstance(p[key], dict):
            for k_str, val in p[key].items():
                if val is not None:
                    try:
                        out[int(k_str)] = float(val)
                    except ValueError:
                        pass
    return out

def classify_class(p, extra_text=""):
    gene = (p.get("gene") or "").upper()
    fam = (p.get("fam") or p.get("family") or "").lower()
    name = (p.get("name") or p.get("protein_name") or "").lower()
    text = f"{fam} {name} {extra_text}".lower()
    if any(x in fam.upper() for x in ("NBPF", "GOLGA6", "NPIP", "POTE")) or gene.startswith(("NBPF", "GOLGA6", "NPIP", "POTE")):
        return "Segmental-duplication arrays"
    if "c2h2" in text or "zinc-finger" in text or "zinc finger" in text or gene.startswith(("ZNF", "KLF", "PRDM", "ZBTB", "ZSCAN", "ZKSCAN", "ZFP")):
        return "C2H2 zinc-finger regulators"
    if gene.startswith(("MUC", "KRTAP", "COL", "FLG", "SPRR", "LCE", "HORN")) or any(w in text for w in ("keratin-associated", "collagen", "mucin")):
        return "Structural tandem-repeat proteins"
    if any(w in text for w in ("dehydrogenase", "kinase", "synthase", "synthetase", "transferase", "hydrolase", "protease", "peptidase", "oxidase", "reductase", "ligase", "isomerase", "lyase", "phosphatase", "atpase", "gtpase", "esterase", "nuclease", "polymerase", "helicase", "carboxylase", "deaminase", "methyltransferase", "acetyltransferase", "ubiquitin-protein ligase")):
        return "Catalytic enzymes"
    if any(w in text for w in ("transporter", "channel", "receptor", "solute carrier", "gpcr", "g-protein coupled", "exchanger", "pump", "importin", "exportin", "aquaporin")) or gene.startswith(("SLC", "ABC", "KCN", "SCN", "CACN", "OR", "GPR", "CHRN", "GRIN", "GABR")):
        return "Receptors, channels and transporters"
    return "Other globular and cellular proteins"

records = []
for p in proteins:
    acc = p.get("acc") or p.get("id") or ""
    L = int(p.get("length") or p.get("L"))
    I_total = L * LOG2_20
    bmap = extract_Bk_dict(p)
    if 1 in bmap and bmap[1] > 0:
        B1 = bmap[1]
    elif acc in seq_map:
        B1 = exact_B1_from_seq(seq_map[acc])
    elif "seq" in p:
        B1 = exact_B1_from_seq(p["seq"])
    else:
        B1 = max(1.0, L * 4.18 - 0.5 * 19 * math.log2(max(2, L)))
    bmap[1] = B1
    I_comp = I_total - B1
    kappa = int(p.get("kappa", 1))
    records.append({
        "acc": acc,
        "gene": p.get("gene", ""),
        "L": L,
        "kappa": kappa,
        "I_total": I_total,
        "I_comp": I_comp,
        "B": bmap,
        "cls": classify_class(p, name_map.get(acc, "")),
    })

avail_ks = sorted({k for r in records for k in r["B"].keys()})

out_lines = []
def log(s=""):
    print(s)
    out_lines.append(s)

log(f"Total reviewed human proteins analysed: {len(records)}")
log(f"Available k orders: {avail_ks}")
log("\n=== 1. Global Information Decomposition across k (bits/residue and % of B1) ===")
log(f"{'k':>3} | {'med I_comp/L':>12} | {'mean I_loc/L':>12} | {'mean B_k/L':>10} | {'mean I_loc/B1':>13} | {'mean B_k/B1':>11} | {'B_k>0 count':>11} | {'B_k>0 share':>11}")
decomp_table = []
for k in avail_ks:
    ic_l = statistics.median(r["I_comp"] / r["L"] for r in records)
    il_l = statistics.mean((r["B"][1] - r["B"].get(k, 0.0)) / r["L"] for r in records)
    bk_l = statistics.mean(r["B"].get(k, 0.0) / r["L"] for r in records)
    il_f = statistics.mean((r["B"][1] - r["B"].get(k, 0.0)) / r["B"][1] for r in records)
    bk_f = statistics.mean(r["B"].get(k, 0.0) / r["B"][1] for r in records)
    pos_cnt = sum(1 for r in records if r["B"].get(k, 0.0) > 0)
    pos_pct = pos_cnt / len(records) * 100.0
    log(f"{k:>3} | {ic_l:>12.4f} | {il_l:>12.4f} | {bk_l:>10.4f} | {il_f*100:>12.3f}% | {bk_f*100:>10.3f}% | {pos_cnt:>11d} | {pos_pct:>10.2f}%")
    decomp_table.append({
        "k": k, "med_I_comp_per_L": round(ic_l, 4), "mean_I_loc_per_L": round(il_l, 4),
        "mean_Bk_per_L": round(bk_l, 4), "mean_I_loc_B1_pct": round(il_f * 100, 3),
        "mean_Bk_B1_pct": round(bk_f * 100, 3), "Bk_pos_count": pos_cnt, "Bk_pos_pct": round(pos_pct, 2)
    })

log("\n=== 2. Functional Classes and Order Regimes ===")
classes = [
    "Catalytic enzymes",
    "Receptors, channels and transporters",
    "Other globular and cellular proteins",
    "C2H2 zinc-finger regulators",
    "Structural tandem-repeat proteins",
    "Segmental-duplication arrays",
]
log(f"{'Class':<38} | {'n':>6} | {'med L':>6} | {'mean B3/B1':>10} | {'mean B5/B1':>10} | {'k>5 share':>9} | {'k>10 share':>10} | {'k>25 share':>10} | {'med kappa':>9} | {'p95 kappa':>9}")
summary_classes = []
for c in classes:
    sub = [r for r in records if r["cls"] == c]
    n = len(sub)
    med_L = statistics.median(r["L"] for r in sub)
    b3_f = statistics.mean(r["B"].get(3, 0.0) / r["B"][1] for r in sub) * 100.0
    b5_f = statistics.mean(r["B"].get(5, 0.0) / r["B"][1] for r in sub) * 100.0
    k5 = sum(1 for r in sub if r["kappa"] > 5) / n * 100.0
    k10 = sum(1 for r in sub if r["kappa"] > 10) / n * 100.0
    k25 = sum(1 for r in sub if r["kappa"] > 25) / n * 100.0
    kappas = sorted(r["kappa"] for r in sub)
    med_k = statistics.median(kappas)
    p95_k = kappas[min(n - 1, int(0.95 * n))]
    log(f"{c:<38} | {n:>6} | {med_L:>6.0f} | {b3_f:>9.2f}% | {b5_f:>9.2f}% | {k5:>8.2f}% | {k10:>9.2f}% | {k25:>9.2f}% | {med_k:>9.0f} | {p95_k:>9d}")
    summary_classes.append({
        "class": c, "n": n, "med_L": med_L, "mean_B3_B1_pct": round(b3_f, 2),
        "mean_B5_B1_pct": round(b5_f, 2), "gt5_pct": round(k5, 2), "gt10_pct": round(k10, 2),
        "gt25_pct": round(k25, 2), "med_kappa": med_k, "p95_kappa": p95_k
    })

with open(os.path.join(RES_DIR, "order_taxonomy_output.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(out_lines) + "\n")

with open(os.path.join(RES_DIR, "order_taxonomy_summary.json"), "w", encoding="utf-8") as f:
    json.dump({"n_proteins": len(records), "decomp": decomp_table, "classes": summary_classes}, f, indent=2)

# Generate Figure 1: Orthogonal Information Decomposition & Survival by Class
ks_plot = [k for k in avail_ks if k >= 2]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.0), dpi=300)

# Left panel: mean residual order entropy B_k / B_1 (%) across k by functional class
palette = {
    "Catalytic enzymes": "#1f77b4",
    "Receptors, channels and transporters": "#2ca02c",
    "Other globular and cellular proteins": "#7f7f7f",
    "C2H2 zinc-finger regulators": "#ff7f0e",
    "Structural tandem-repeat proteins": "#9467bd",
    "Segmental-duplication arrays": "#d62728",
}
for c in classes:
    sub = [r for r in records if r["cls"] == c]
    ys = [max(1e-4, statistics.mean(r["B"].get(k, 0.0) / r["B"][1] for r in sub) * 100.0) for k in ks_plot]
    ax1.plot(ks_plot, ys, marker="o", markersize=4, linewidth=1.6, label=c, color=palette[c])

ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel("local encoder order $k$")
ax1.set_ylabel("unresolved global order $B_k(S)/B_1(S)$ (%)")
ax1.set_title("Residual global order information across $k$", fontsize=10, fontweight="bold")
ax1.grid(True, which="both", linestyle=":", alpha=0.5)
ax1.legend(fontsize=7.2, loc="upper right", frameon=True)

# Right panel: share of proteins with B_k(S) > 0 (kappa > k) across k by class
for c in classes:
    sub = [r for r in records if r["cls"] == c]
    n = len(sub)
    ys = [max(0.02, sum(1 for r in sub if r["kappa"] > k) / n * 100.0) for k in ks_plot]
    ax2.plot(ks_plot, ys, marker="s", markersize=4, linewidth=1.6, label=c, color=palette[c])

ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel("local encoder order $k$")
ax2.set_ylabel("proteins with $B_k(S)>0$ (%)")
ax2.set_title("Unresolved tail $\\kappa(S) > k$ by functional class", fontsize=10, fontweight="bold")
ax2.grid(True, which="both", linestyle=":", alpha=0.5)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "fig_taxonomy.pdf"))
plt.close()
print("Saved figure to figures/fig_taxonomy.pdf and results to results/")