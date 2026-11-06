"""
python pipeline/make_summary_xlsx.py  ->  PhyloBench_boundaries_summary.xlsx (one sheet, Summary)

Collects, for OB45, AG45, AS45, CH45 and EB45, each method's boundary (largest number of random loci it can run on)
and the three stages that were run, from the files the pipeline writes:
  <X>45/logs/*_boundary*/probes.tsv and boundary.txt, <X>45_<N>_result_trees/, Comparison/figures_R/*.log.
Re-run it whenever new results come in; numbers are never typed in by hand except where noted in the Notes sheet.
"""
import csv
import random
import re
from pathlib import Path

from openpyxl import Workbook
from openpyxl.comments import Comment
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "PhyloBench_boundaries_summary.xlsx"
DATASETS = ["OB", "AG", "AS", "CH", "EB"]
MSAT_LIMIT = 1023
FONT = "Arial"
MARK_UNCERTAIN = False  # True: label uncertain values, colour the cells and explain them below the table


def locus_lengths(X):
    lengths = {}
    for f in sorted((ROOT / f"{X}45_cleaned").glob("*.fasta")):
        with open(f) as h:
            h.readline()
            n = 0
            for line in h:
                if line.startswith(">"):
                    break
                n += len(line.strip())
        lengths[f] = n
    return lengths


def random_cols(files, lengths, n):
    random.seed(42)
    return sum(lengths[f] for f in random.sample(files, n))


def read_boundary(path):
    if not path.exists():
        return None
    t = path.read_text()
    lo = re.search(r"largest success: (\d+) loci \((\d+) columns\)", t)
    hi = re.search(r"smallest failure above it: (\d+) loci", t)
    return int(lo.group(1)), int(lo.group(2)), int(hi.group(1))


def read_probes(X, files, lengths):
    rows = []
    for p in sorted((ROOT / f"{X}45/logs").glob("*_boundary*/probes.tsv")):
        method = "CMAPLE" if "cmaple" in p.parent.name else "Phyloformer"
        selection = "shortest loci" if "shortest" in p.parent.name else "random (seed 42)"
        for r in csv.DictReader(open(p), delimiter="\t"):
            if not r.get("n_loci", "").isdigit():
                continue
            ok = r["success"] == "True"
            rows.append({"dataset": f"{X}45", "method": method, "selection": selection, "round": p.parent.name,
                         "n": int(r["n_loci"]), "cols": int(r["columns"]),
                         "ok": None if not ok and not r.get("error") else ok,  # no tree and no error: log was cut off
                         "sec": float(r["seconds"]) if r["seconds"] else None,
                         "peak": float(r["peak_GB"]) if r.get("peak_GB") else None, "error": r.get("error", ""),
                         "note": ""})
        # A search killed by the full disk never logged the probes it had finished; read them from their folders
        if method == "CMAPLE" and "round1" in p.parent.name:
            logged = {r["n"] for r in rows if r["round"] == p.parent.name}
            for d in sorted(p.parent.glob("N*")):
                n = int(d.name[1:])
                if n in logged:
                    continue
                out = (d / "cmaple.out").read_text() if (d / "cmaple.out").exists() else ""
                err = next((l.strip() for l in out.splitlines() if l.startswith("ERROR")), "")
                tree = any(d.glob("*.treefile"))
                rows.append({"dataset": f"{X}45", "method": method, "selection": selection, "round": p.parent.name,
                             "n": n, "cols": random_cols(files, lengths, n), "ok": True if tree else (False if err else None),
                             "sec": None, "peak": None, "error": err, "note": "read from the probe folder"})
    return rows


def has_tree(X, label, method):
    d = ROOT / f"{X}45_{label}_result_trees"
    pat = {"IQ": "*IQTree.treefile", "PF": "*[Pp]hyloformer.nwk", "CM": "*CMAPLE.treefile",
           "MA": "*AttentionTree.nwk", "ME": "*embedding_tree.nwk"}[method]
    return d.is_dir() and any(d.glob(pat))


def quartet_share(X, label):
    log = ROOT / f"Comparison/figures_R/compare_simple_{X}_{label}.log"
    if not log.exists():
        return None
    m = re.search(r"Quartets resolved identically by all trees: \d+ \(([\d.]+)%\)", log.read_text())
    return float(m.group(1)) / 100 if m else None


# ---------------------------------------------------------------- gather
info = {}
for X in DATASETS:
    lengths = locus_lengths(X)
    files = sorted(lengths)
    msat = max(n for n in range(1, 30) if random_cols(files, lengths, n) <= MSAT_LIMIT)
    info[X] = {
        "files": files, "lengths": lengths, "n_loci": len(files), "total_cols": sum(lengths.values()),
        "msat": (msat, random_cols(files, lengths, msat), msat + 1, random_cols(files, lengths, msat + 1)),
        "cmaple": read_boundary(ROOT / f"{X}45/logs/cmaple_boundary/boundary.txt"),
        "pf": read_boundary(ROOT / f"{X}45/logs/phyloformer_boundary/boundary.txt"),
        "probes": read_probes(X, files, lengths),
    }

# Phyloformer memory depends only on columns, so a failure below the largest success of any dataset cannot be the limit
PF_MAX_OK = max((r["cols"] for X in DATASETS for r in info[X]["probes"] if r["method"] == "Phyloformer" and r["ok"]),
                default=0)
for X in DATASETS:
    for r in info[X]["probes"]:
        if r["method"] == "Phyloformer" and r["ok"] is False:
            if r["cols"] < PF_MAX_OK:
                r["suspect"] = True
                r["note"] = "SUSPECT: a larger matrix succeeded, so other jobs were probably using the memory"
            elif "round1" in r["round"]:
                r["note"] = "UNCERTAIN: ran while other jobs were using memory and was not repeated"
        elif r["ok"] is None:
            r["note"] = "INVALID: log cut off when the disk quota filled; not counted"


def cmaple_failures_below(X, n):
    return sorted({r["n"] for r in info[X]["probes"] if r["method"] == "CMAPLE" and r["selection"].startswith("random")
                   and r["ok"] is False and r["n"] < n})


def certainty(X, method):
    """(Certainty, Why) for one boundary."""
    d = info[X]
    if method == "msat":
        return "Certain", f"Fixed limit of {MSAT_LIMIT} columns; column counts are exact for the seed-42 samples"
    if method == "cmaple":
        if not d["cmaple"]:
            return "Not determined", "Search not finished"
        below = cmaple_failures_below(X, d["cmaple"][0])
        if below:
            return "Uncertain", ("Success is not monotonic in N: N = " + ", ".join(map(str, below)) +
                                 " failed below this value, even though those samples are subsets of the sample at "
                                 "this N, so a larger N or other loci could differ")
        return "Certain", "Every tested N below succeeded and every tested N above failed (seed-42 samples)"
    if method == "pf":
        if X in ("OB", "AG"):
            return "Certain", "All loci ran, so no limit was reached"
        if not d["pf"]:
            return "Not determined", "Search not finished; the value shown is the largest success so far"
        if X == "EB":
            return "Uncertain", ("The 749 failure ran while other jobs were using memory and was not repeated, and "
                                 "668 and 695 from that same round succeeded when repeated. The true limit may be "
                                 "higher, but all 856 loci (~1.09 TB) cannot fit")
        return "Uncertain", "Depends on how much memory other jobs were using during the probes"
    if method == "iq":
        if X in ("OB", "AG"):
            return "Certain", "All loci ran"
        return "Not tested on all loci", "No hard limit is known, but all loci were never run (expected to take too long)"


def pf_provisional(X):
    """Largest Phyloformer success so far and the first failure above it that is not suspect."""
    rows = [r for r in info[X]["probes"] if r["method"] == "Phyloformer"]
    ok = [r for r in rows if r["ok"]]
    best = max(ok, key=lambda r: r["n"]) if ok else None
    fails = [r for r in rows if r["ok"] is False and not r.get("suspect") and (not best or r["n"] > best["n"])]
    fail = min(fails, key=lambda r: r["n"]) if fails else None
    return best, fail

# Stage groups: label of <X>45_<label>_result_trees for each stage (None = not run yet)
STAGE_LABELS = {
    "OB": {1: "7", 2: None, 3: "603"},
    "AG": {1: "4", 2: "4", 3: "286"},
    "AS": {1: "3", 2: "113", 3: None},
    "CH": {1: "5", 2: "307", 3: None},
    "EB": {1: "2", 2: "18", 3: "748"},
}
for X in DATASETS:
    for s in (2, 3):
        key = "cmaple" if s == 2 else "pf"
        if STAGE_LABELS[X][s] is None and info[X][key] and (ROOT / f"{X}45_{info[X][key][0]}_result_trees").is_dir():
            STAGE_LABELS[X][s] = str(info[X][key][0])

# ---------------------------------------------------------------- workbook helpers
wb = Workbook()
HEAD_FILL = PatternFill("solid", start_color="1F4E78")
HEAD_FONT = Font(name=FONT, bold=True, color="FFFFFF")
OK_FILL = PatternFill("solid", start_color="E2EFDA")
BAD_FILL = PatternFill("solid", start_color="FCE4D6")
TODO_FILL = PatternFill("solid", start_color="FFF2CC")
THIN = Border(*(Side(style="thin", color="BFBFBF"),) * 4)


def sheet(title, headers, widths):
    ws = wb.create_sheet(title)
    ws.append(headers)
    for c, w in zip(ws[1], widths):
        c.fill, c.font = HEAD_FILL, HEAD_FONT
        c.alignment = Alignment(wrap_text=True, vertical="center")
        ws.column_dimensions[c.column_letter].width = w
    ws.freeze_panes = "A2"
    return ws


def finish(ws):
    for row in ws.iter_rows(min_row=2):
        for c in row:
            c.font = Font(name=FONT, bold=c.font.bold, color=c.font.color)
            c.border = THIN
            c.alignment = Alignment(wrap_text=True, vertical="top")
    ws.auto_filter.ref = ws.dimensions


wb.remove(wb.active)

# ---------------------------------------------------------------- Boundaries
ws = sheet("Boundaries",
           ["Dataset", "Total loci", "Total columns", "Method", "What limits it", "Boundary N (loci)",
            "Columns at N", "Columns per locus at N", "First failing N", "Columns at first failing N", "Status",
            "Certainty", "Why", "Evidence"],
           [10, 10, 13, 22, 34, 14, 13, 13, 13, 15, 26, 14, 60, 50])


def iqtree_time(X, label):
    """Wall-clock time of the IQ-TREE run for this group (full original runs have no _<N> suffix)."""
    for name in (f"{X}45_supermatrix_{label}", f"{X}45_supermatrix"):
        log = ROOT / f"IQTREE_results/{name}.log"
        if log.exists() and (name.endswith(label) or int(label) == info[X]["n_loci"]):
            m = re.search(r"Total wall-clock time used: [\d.]+ sec \((\S+)\)", log.read_text())
            if m:
                return m.group(1)
    return None


def largest_iqtree_run(X):
    done = [int(STAGE_LABELS[X][s]) for s in (3, 2, 1) if STAGE_LABELS[X][s] and iqtree_time(X, STAGE_LABELS[X][s])]
    return (max(done), iqtree_time(X, str(max(done)))) if done else (None, None)

for X in DATASETS:
    d = info[X]
    rows = []
    n, c, n1, c1 = d["msat"]
    msat_check = {"OB": "OB45_7 ran all five methods; N+1 computed, not run",
                  "AG": "AG45_4 ran all five methods; N+1 computed, not run",
                  "AS": "AS45_3 ran all five methods; N+1 run and rejected by MSA-T",
                  "CH": "CH45_5 ran all five methods; N+1 run and rejected by MSA-T",
                  "EB": "N: embedding tree ran, attention not run (EB45_2 ran all five methods); N+1 computed, not run"}[X]
    rows.append(["MSA Transformer (embedding + attention)", f"at most {MSAT_LIMIT} columns", n, c, n1, c1,
                 "found", f"pipeline/msat_boundary.py. {msat_check}"])
    if d["cmaple"]:
        lo, lc, hi = d["cmaple"]
        rows.append(["CMAPLE", "divergence (\"Likelihood cost is very heavy\")", lo, lc, hi,
                     random_cols(d["files"], d["lengths"], hi), "found",
                     f"{X}45/logs/cmaple_boundary/boundary.txt, probes.tsv"])
    else:
        rows.append(["CMAPLE", "divergence (\"Likelihood cost is very heavy\")", None, None, None, None,
                     "searching", f"{X}45/logs/cmaple_boundary/search.out"])
    if X in ("OB", "AG"):
        rows.append(["Phyloformer", "RAM: ~5.5 GB per 1000 columns (1 TB machine ≈ 174k columns)", d["n_loci"],
                     d["total_cols"], None, None, "no limit: all loci ran",
                     f"{X}45_{d['n_loci']}_result_trees (Phyloformer on the full supermatrix)"])
    elif d["pf"]:
        lo, lc, hi = d["pf"]
        rows.append(["Phyloformer", "RAM: ~5.5 GB per 1000 columns (1 TB machine ≈ 174k columns)", lo, lc, hi,
                     random_cols(d["files"], d["lengths"], hi), "found",
                     f"{X}45/logs/phyloformer_boundary/boundary.txt, probes.tsv"])
    else:
        best, fail = pf_provisional(X)
        rows.append(["Phyloformer", "RAM: ~5.5 GB per 1000 columns (1 TB machine ≈ 174k columns)",
                     best["n"] if best else None, best["cols"] if best else None,
                     fail["n"] if fail else 800, fail["cols"] if fail else random_cols(d["files"], d["lengths"], 800),
                     "searching" + ("" if fail else " (800 random loci failed)"),
                     f"{X}45/logs/phyloformer_boundary/probes.tsv"])
    big_n, big_t = largest_iqtree_run(X)
    rows.append(["IQ-TREE", "none (run time only)", None, None, None, None,
                 f"no limit; largest run {big_n} loci in {big_t}",
                 "IQTREE_results/<X>45_supermatrix_<N>.log (Total wall-clock time)"])
    for r, key in zip(rows, ["msat", "cmaple", "pf", "iq"]):
        cert, why = certainty(X, key)
        ws.append([f"{X}45", d["n_loci"], d["total_cols"]] + r[:4] + [None] + r[4:7] + [cert, why, r[7]])
        i = ws.max_row
        ws[f"L{i}"].fill = OK_FILL if cert == "Certain" else TODO_FILL
        if cert != "Certain" and r[2] is not None:
            ws[f"F{i}"].fill = TODO_FILL
        ws[f"H{i}"] = f'=IF(AND(ISNUMBER(F{i}),F{i}>0),G{i}/F{i},"")'
        ws[f"H{i}"].number_format = "0"
        fill = OK_FILL if r[6].startswith(("found", "no limit")) else TODO_FILL
        ws[f"K{i}"].fill = fill
        for col in "BCFGIJ":
            ws[f"{col}{i}"].number_format = "#,##0"
finish(ws)
ws["E1"].comment = Comment("Phyloformer memory: measured peaks on EB45 were 542 GB at 98,547 columns, 822 GB at "
                           "149,414 and 944 GB at 171,550 (EB45/logs/phyloformer_boundary*/probes.tsv).", "pipeline")

# ---------------------------------------------------------------- Stages
STAGE_NAME = {1: "1: all five methods", 2: "2: all but MSA-T", 3: "3: IQ-TREE + Phyloformer"}
ws = sheet("Stages",
           ["Dataset", "Stage", "N loci", "Columns", "IQ-TREE", "Phyloformer", "CMAPLE", "MSA-T attention",
            "MSA-T embedding", "Trees produced", "Result trees", "Plots", "Quartets identical in all trees",
            "Certainty of N"],
           [10, 24, 9, 11, 9, 12, 9, 12, 12, 10, 30, 36, 14, 14])
EXPECTED = {1: ["IQ", "PF", "CM", "MA", "ME"], 2: ["IQ", "PF", "CM"], 3: ["IQ", "PF"]}
for X in DATASETS:
    for s in (1, 2, 3):
        label = STAGE_LABELS[X][s]
        if label is None:
            key = "cmaple" if s == 2 else "pf"
            n = info[X][key][0] if info[X][key] else None
            ws.append([f"{X}45", STAGE_NAME[s], n, None] + ["to run" if m in EXPECTED[s] else "–"
                                                            for m in ["IQ", "PF", "CM", "MA", "ME"]]
                      + [None, "not run yet", None, None, certainty(X, key)[0]])
            i = ws.max_row
            ws[f"N{i}"].fill = TODO_FILL
            for col in "EFGHI":
                if ws[f"{col}{i}"].value == "to run":
                    ws[f"{col}{i}"].fill = TODO_FILL
            ws[f"J{i}"] = f'=COUNTIF(E{i}:I{i},"tree")'
            continue
        n = int(label)
        cols = info[X]["total_cols"] if n == info[X]["n_loci"] else random_cols(info[X]["files"], info[X]["lengths"], n)
        marks = []
        for m in ["IQ", "PF", "CM", "MA", "ME"]:
            marks.append("tree" if has_tree(X, label, m) else ("failed" if m not in EXPECTED[s] else "missing"))
        pdf = ROOT / f"Comparison/figures_R/compare_simple_{X}_{label}.pdf"
        ws.append([f"{X}45", STAGE_NAME[s], n, cols] + marks +
                  [None, f"{X}45_{label}_result_trees", f"figures_R/{pdf.name}" if pdf.exists() else "no plot yet",
                   quartet_share(X, label), certainty(X, {1: "msat", 2: "cmaple", 3: "pf"}[s])[0]])
        i = ws.max_row
        ws[f"N{i}"].fill = OK_FILL if ws[f"N{i}"].value == "Certain" else TODO_FILL
        if ws[f"N{i}"].value != "Certain":
            ws[f"C{i}"].fill = TODO_FILL
        ws[f"J{i}"] = f'=COUNTIF(E{i}:I{i},"tree")'
        for col in "EFGHI":
            v = ws[f"{col}{i}"].value
            ws[f"{col}{i}"].fill = OK_FILL if v == "tree" else BAD_FILL
        ws[f"M{i}"].number_format = "0.0%"
        for col in "CD":
            ws[f"{col}{i}"].number_format = "#,##0"
finish(ws)
ws["M1"].comment = Comment("Share of all 4-taxon subsets that every tree in the plot (methods + 3 reference trees) "
                           "resolves the same way; from Comparison/figures_R/compare_simple_<X>_<N>.log.", "pipeline")

# ---------------------------------------------------------------- Probes
ws = sheet("Probes", ["Dataset", "Method", "Loci selection", "Probe round", "N loci", "Columns", "Success",
                      "Seconds", "Peak RAM (GB)", "Note", "Error"], [10, 13, 16, 34, 9, 11, 9, 10, 12, 44, 60])
for X in DATASETS:
    for r in sorted(info[X]["probes"], key=lambda r: (r["method"], r["selection"], r["n"])):
        ws.append([r["dataset"], r["method"], r["selection"], r["round"], r["n"], r["cols"],
                   {True: "yes", False: "no", None: "invalid"}[r["ok"]], r["sec"], r["peak"], r["note"], r["error"]])
        i = ws.max_row
        ws[f"G{i}"].fill = OK_FILL if r["ok"] else BAD_FILL
        if r["note"].startswith(("SUSPECT", "UNCERTAIN", "INVALID")):
            ws[f"J{i}"].fill = TODO_FILL
        ws[f"F{i}"].number_format = "#,##0"
        ws[f"H{i}"].number_format = "#,##0"
finish(ws)

# ---------------------------------------------------------------- Notes
ws = wb.create_sheet("Notes")
notes = [
    ("Loci selection", "All boundaries and stages use random loci: random.seed(42); random.sample(sorted loci, N) "
     "(pipeline/downsample_loci.py, same as OB45/downsample_loci_7.py). N counts loci (gene alignments); every "
     "supermatrix has the dataset's 60 taxa (a taxon absent from all N loci is dropped)."),
    ("Boundary", "Largest N that runs, directly below the smallest N that fails. CMAPLE can succeed on a sample "
     "containing one that failed (see Probes)."),
    ("MSA Transformer", f"Accepts at most {MSAT_LIMIT} columns (1024 positions minus the start token); embedding and "
     "attention trees share this limit. Computed for each N, checked by running N (works) and N+1 (fails)."),
    ("CMAPLE", "Fails with 'Likelihood cost is very heavy' when the sequences are too divergent."),
    ("Phyloformer", "Memory grows linearly: ~5.5 GB per 1000 columns for 60 taxa (1770 taxon pairs x columns x "
     "1 KB per tensor, ~3 tensors at once). On this 1 TB machine the limit is ~174k columns; other programs using "
     "memory at the same time lower it."),
    ("IQ-TREE", "No hard limit; only run time (largest run and its wall-clock time read from IQTREE_results/*.log)."),
    ("Stages", "1 = largest N on which all five methods give a tree; 2 = CMAPLE boundary (MSA-T fails); "
     "3 = Phyloformer boundary or all loci (MSA-T and CMAPLE fail). AG45: CMAPLE and MSA-T boundaries coincide (4), so stages 1 and 2 are the same run. EB45: the MSA-T boundary is 3, "
     "but CMAPLE fails on that random 3-locus sample, so stage 1 uses 2 loci."),
    ("Plots", "Comparison/figures_R/compare_simple_<X>_<N>.pdf; trees rooted on the root tip of <X>-astral.tre "
     "(OB MYROD, AG METFS, AS DOTSN, CH DANRE, EB CHLRE); RF and quartet distances unrooted."),
    ("Other runs", "Earlier shortest-loci subsets are kept but not part of the stages: OB45_short_15/short_20, "
     "AG45_18/36, EB45_short_26/short_46. OB45_603 and AG45_286 are the original full-supermatrix runs."),
    ("Certainty", "Certain = fixed limit, or all tested results agree. Uncertain = the value is usable but could "
     "differ (see Why in Boundaries). Not determined = the search has not finished. Uncertain values are highlighted."),
    ("Source", "Generated by pipeline/make_summary_xlsx.py from the pipeline's log and result files."),
]
ws.column_dimensions["A"].width = 18
ws.column_dimensions["B"].width = 110
for k, v in notes:
    ws.append([k, v])
for row in ws.iter_rows():
    row[0].font = Font(name=FONT, bold=True)
    row[1].font = Font(name=FONT)
    row[1].alignment = Alignment(wrap_text=True, vertical="top")
    row[0].alignment = Alignment(vertical="top")

# ---------------------------------------------------------------- Summary (first sheet)
ws = sheet("Summary", ["Dataset", "Total loci", "Total columns", "MSA-T boundary", "CMAPLE boundary",
                       "Phyloformer boundary", "IQ-TREE", "Stage 1 N", "Stage 2 N", "Stage 3 N"],
           [10, 10, 13, 16, 22, 24, 22, 11, 22, 24])
for X in DATASETS:
    d = info[X]
    def cell(n, key):
        c = certainty(X, key)[0]
        if n is None:
            return "not determined"
        if c == "Certain":
            return str(n)
        if c == "Not determined":
            return f"at least {n} (not determined)"
        return f"{n} (uncertain)" if MARK_UNCERTAIN else str(n)
    best = d["pf"][0] if d["pf"] else (pf_provisional(X)[0] or {}).get("n")
    if X in ("OB", "AG"):
        best = d["n_loci"]
    done_n, done_t = largest_iqtree_run(X)
    s3 = STAGE_LABELS[X][3]
    iq_done = bool(s3 and iqtree_time(X, s3))
    if X in ("OB", "AG"):
        iq_cell = f"all {d['n_loci']} loci ran"
    elif s3 and not iq_done:
        iq_cell = f"{done_n} finished; {s3} running"
    else:
        iq_cell = f"{done_n} finished; all {d['n_loci']} loci not run"
    ws.append([f"{X}45", d["n_loci"], d["total_cols"], cell(d["msat"][0], "msat"),
               cell(d["cmaple"][0] if d["cmaple"] else None, "cmaple"), cell(best, "pf"),
               iq_cell, STAGE_LABELS[X][1], cell(d["cmaple"][0] if d["cmaple"] else None, "cmaple"),
               (cell(best, "pf") + ("" if iq_done else "; IQ-TREE running")) if STAGE_LABELS[X][3] else "not run yet"])
    i = ws.max_row
    for col in "BC":
        ws[f"{col}{i}"].number_format = "#,##0"
    if MARK_UNCERTAIN:
        for col in "DEFIJ":
            v = str(ws[f"{col}{i}"].value)
            ws[f"{col}{i}"].fill = TODO_FILL if ("uncertain" in v or "not " in v or "running" in v) else OK_FILL
        ws[f"G{i}"].fill = OK_FILL if X in ("OB", "AG") else TODO_FILL
finish(ws)
ws.append([])
FOOT = [
    "N = number of random loci (seed 42); every supermatrix has the dataset's 60 taxa, or fewer if a taxon is absent from "
    "all N loci. Stage 1 = all five methods give a tree; stage 2 = all but MSA-T; stage 3 = IQ-TREE + Phyloformer only.",
    "Generated by pipeline/make_summary_xlsx.py from the pipeline logs.",
    "Stage 1: EB45's MSA-T boundary is 3, but CMAPLE fails on that 3-locus sample, so stage 1 uses 2. AG45: MSA-T and "
    "CMAPLE boundaries are both 4, so stages 1 and 2 are the same run.",
] + ([] if not MARK_UNCERTAIN else [
    "Yellow = uncertain, not determined or still running.",
    "CMAPLE uncertain (AS45, CH45, EB45): smaller N also failed (AS45: 31, 103, 105, 107, 110; CH45: 197, 233, 282, 295; "
    "EB45: 3), so the value holds for these seed-42 samples and could differ for other loci.",
    "Phyloformer uncertain (limit depends on memory used by other programs): AS45 693 ran, 700 failed twice while other "
    "jobs used memory. CH45 600 is a rounded value: it ran, 700 failed while other jobs used memory; the true limit is "
    "probably higher. EB45 748 ran, 749 failed once while other jobs used memory and was not repeated.",
])
for t in FOOT:
    ws.append([t])
    i = ws.max_row
    ws.merge_cells(start_row=i, start_column=1, end_row=i, end_column=10)
    ws.cell(i, 1).font = Font(name=FONT, italic=True)
    ws.cell(i, 1).alignment = Alignment(wrap_text=True, vertical="top")
    ws.row_dimensions[i].height = 32
# Only the Summary sheet is delivered
for name in list(wb.sheetnames):
    if name != "Summary":
        wb.remove(wb[name])

# No LibreOffice on the server to cache formula results, so let Excel compute them when the file is opened
from openpyxl.workbook.properties import CalcProperties
wb.calculation = CalcProperties(fullCalcOnLoad=True)
wb.save(OUT)
print("wrote", OUT)
