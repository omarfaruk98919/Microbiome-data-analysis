#!/usr/bin/env python3
"""
ARG Processing Pipeline
=======================
Processes antibiotic resistance gene (ARG) data through the following steps:
  1. Deduplicate hits from a multi-database .tab file (position-based overlap)
  2. Merge duplicate GENE entries across databases
  3. Fill missing RESISTANCE annotations using gene-level lookup
  4. Apply manual gene remappings and direct value injections
  5. Map resistance labels to standardised drug classes
  6. Classify entries (single-class / dual-class / multidrug) and explode rows

Usage
-----
    python process_ARGs.py

Inputs / Outputs
----------------
    Input  : alldb_rare_results.tab              (raw ABRicate-style tab file)
    Output : deduplicated_ARGs_rare.csv           (step 1)
             deduplicated_ARGs_rare_casemerged.csv (step 2)
             deduplicated_ARGs_rare_resistance_filled.csv (step 3)
             deduplicated_ARGs_rare_resistance_filled_final.csv (step 4)
             deduplicated_ARGs_rare_mapped.csv    (step 5 + 6, final output)
"""

import csv
import re
import pandas as pd

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

INPUT_TAB   = "alldb_rare_results.tab"
OVERLAP_TOL = 100                          # nt tolerance for positional overlap
DB_PRIORITY = ["CARD", "NCBI", "RESFINDER", "ARGANNOT", "MEGARES"]
PRIORITY_MAP = {db: i for i, db in enumerate(DB_PRIORITY)}

# Resistance label → standardised drug class
DICT_MAP = {
    # Aminoglycoside
    "aminoglycoside": "Aminoglycoside",
    "streptomycin":   "Aminoglycoside",
    "gentamicin":     "Aminoglycoside",
    "kanamycin":      "Aminoglycoside",
    "amikacin":       "Aminoglycoside",
    "tobramycin":     "Aminoglycoside",
    "neomycin":       "Aminoglycoside",
    "plazomicin":     "Aminoglycoside",

    # Beta-lactam & subclasses
    "beta lactam":          "Beta-lactam",
    "betalactam":           "Beta-lactam",
    "beta-lactam":          "Beta-lactam",
    "penicillin":           "Beta-lactam",
    "penicillin beta lactam":  "Beta-lactam",
    "penicillin beta-lactam":  "Beta-lactam",
    "cephalosporin":        "Beta-lactam",
    "carbapenem":           "Beta-lactam",
    "monobactam":           "Beta-lactam",
    "oxacephem":            "Beta-lactam",
    "penem":                "Beta-lactam",
    "penam":                "Beta-lactam",
    "cephamycin":           "Beta-lactam",

    # Fluoroquinolone
    "fluoroquinolone": "Fluoroquinolone",
    "ciprofloxacin":   "Fluoroquinolone",
    "levofloxacin":    "Fluoroquinolone",
    "norfloxacin":     "Fluoroquinolone",

    # Fosfomycin
    "fosfomycin":      "Fosfomycin",
    "phosphonic acid": "Fosfomycin",

    # MLSB
    "macrolide":     "MLSB",
    "lincosamide":   "MLSB",
    "streptogramin": "MLSB",
    "streptogramin b": "MLSB",
    "clindamycin":   "MLSB",
    "erythromycin":  "MLSB",
    "azithromycin":  "MLSB",
    "tylosin":       "MLSB",
    "spiramycin":    "MLSB",
    "telithromycin": "MLSB",

    # Phenicol
    "phenicol":      "Phenicol",
    "chloramphenicol": "Phenicol",
    "florfenicol":   "Phenicol",

    # Rifamycin
    "rifamycin":  "Rifamycin",
    "rifampin":   "Rifamycin",
    "rifampicin": "Rifamycin",

    # Sulfonamide
    "sulfonamide":  "Sulfonamide",
    "sulphonamide": "Sulfonamide",

    # Tetracycline (incl. glycylcyclines)
    "tetracycline":  "Tetracycline",
    "doxycycline":   "Tetracycline",
    "minocycline":   "Tetracycline",
    "glycylcycline": "Tetracycline",
    "tigecycline":   "Tetracycline",

    # Peptide
    "peptide":   "Peptide",
    "bacitracin": "Peptide",
    "polymyxin":  "Peptide",

    # Nucleoside
    "nucleoside": "Nucleoside",

    # Trimethoprim
    "diaminopyrimidine": "Trimethoprim",
    "trimethoprim":      "Trimethoprim",

    # Glycopeptide
    "glycopeptide": "Glycopeptide",
    "bleomycin":    "Glycopeptide",

    # Aminocoumarin
    "aminocoumarin": "Aminocoumarin",

    # Oxazolidinone
    "oxazolidinone": "Oxazolidinone",
    "oxazolidin":    "Oxazolidinone",

    # Nitrofuran
    "nitofuran":     "Nitrofuran",
    "nitroimidazole": "Nitrofuran",

    # Pleuromutilin
    "pleuromutilin": "Pleuromutilin",

    # Biocides → drop
    "disinfecting agents and antiseptics": "_DROP_BIOCIDE",
}

# Individual items to silently drop
DROP_ITEMS = {
    "triclosan",
    "acridine_dye",
    "benzalconium_chloride",
    "benzalkonium_chloride",
    "rhodamine",
    "fucidic_acid",
}

# Manual gene-to-gene resistance remapping  (target gene copies source gene's resistance)
MANUAL_MAPPING = {
    "CPXAR": "CPXA",
    "NIME_1": "NIME",
}

# Direct resistance value injection for specific genes
MANUAL_VALUES = {
    "EMRD": "Phenicol",
}


# ─────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def db_rank(db):
    return PRIORITY_MAP.get(str(db).strip().upper(), 10**9)


def split_dbs(s):
    """Split a semicolon- or comma-separated DB string into a list of tokens."""
    return [t.strip().upper() for t in str(s).replace(",", ";").split(";") if t.strip()]


def best_db_rank(cell):
    toks = split_dbs(cell)
    return min(db_rank(t) for t in toks) if toks else 10**9


def cluster_mask(row, d):
    """Return boolean mask for rows that overlap with `row` within OVERLAP_TOL."""
    base = (
        (d["FILE"] == row["FILE"]) &
        (d["SEQUENCE"] == row["SEQUENCE"]) &
        (d["STRAND"] == row["STRAND"])
    )
    return base & (
        ((d["START"] - row["START"]).abs() <= OVERLAP_TOL) |
        ((d["END"]   - row["END"]).abs()   <= OVERLAP_TOL)
    )


def merge_db_tokens(db_series):
    """Union all unique DB tokens across a series of semicolon-separated strings."""
    tokens = set()
    for entry in db_series.dropna():
        for tok in str(entry).split(";"):
            t = tok.strip()
            if t:
                tokens.add(t)
    return ";".join(sorted(tokens, key=lambda x: (db_rank(x), x)))


def merge_resistance_tokens(res_series):
    """Union all unique resistance tokens across a series."""
    all_res = set()
    for entry in res_series.dropna():
        if str(entry).lower() in {"none", "nan", ""}:
            continue
        for tok in str(entry).split(";"):
            t = tok.strip()
            if t:
                all_res.add(t)
    return ";".join(sorted(all_res)) if all_res else None


def map_resistance_list(text):
    """Map a semicolon-separated resistance string through DICT_MAP, dropping biocides."""
    if pd.isna(text):
        return text
    mapped = []
    for item in [x.strip().lower() for x in str(text).split(";")]:
        if item in DROP_ITEMS:
            continue
        if DICT_MAP.get(item) == "_DROP_BIOCIDE":
            continue
        mapped.append(DICT_MAP.get(item, item))
    return ";".join(sorted(set(mapped)))


def classify_resistance(mapped_text):
    """
    Return (classification_label, [classes]) where:
      - classification_label is 'multidrug' (>=3), None (1-2), or '' (empty)
      - classes is the sorted unique list of drug classes present
    """
    if pd.isna(mapped_text) or mapped_text.strip() == "":
        return "", []
    classes = sorted({x.strip() for x in mapped_text.split(";") if x.strip()})
    if len(classes) >= 3:
        return "multidrug", classes
    return None, classes


# ─────────────────────────────────────────────────────────────────────────────
# STEP 1 — Deduplicate by positional overlap
# ─────────────────────────────────────────────────────────────────────────────

def step1_deduplicate(input_tab):
    print("\n── Step 1: Positional deduplication ──")
    df = pd.read_csv(input_tab, sep="\t", dtype=str)

    for c in ["START", "END"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    for c in ["FILE", "SEQUENCE", "STRAND", "DATABASE"]:
        if c in df.columns:
            df[c] = df[c].astype(str).str.strip()
    df["DATABASE"] = df["DATABASE"].str.upper()

    df = df.sort_values(
        by=[c for c in ["FILE", "SEQUENCE", "START"] if c in df.columns],
        kind="mergesort"
    ).reset_index(drop=True)

    used = pd.Series(False, index=df.index)
    merged_rows = []

    for i in df.index:
        if used[i]:
            continue
        sub = df.loc[cluster_mask(df.loc[i], df) & (~used)].copy()

        merged_db = merge_db_tokens(sub["DATABASE"])

        sub["BEST_DB_RANK"]  = sub["DATABASE"].map(best_db_rank)
        sub["NUM_DB_TOKENS"] = sub["DATABASE"].map(lambda x: len(split_dbs(x)))
        rep = sub.sort_values(
            by=["BEST_DB_RANK", "NUM_DB_TOKENS", "START"],
            ascending=[True, True, True],
            kind="mergesort"
        ).iloc[0].copy()
        rep["DATABASE"] = merged_db

        merged_rows.append(rep)
        used.loc[sub.index] = True

    final_df = pd.DataFrame(merged_rows).sort_values(
        by=[c for c in ["FILE", "SEQUENCE", "START"] if c in df.columns],
        kind="mergesort"
    ).reset_index(drop=True)

    # Protect GAPS from Excel date-parsing
    if "GAPS" in final_df.columns:
        final_df["GAPS"] = final_df["GAPS"].astype(str).str.strip()
        final_df["GAPS"] = final_df["GAPS"].map(
            lambda s: s if s.startswith("*") or not re.fullmatch(r"\d{1,4}/\d{1,4}", s) else f"*{s}"
        )

    out = "deduplicated_ARGs_rare.csv"
    final_df.to_csv(out, index=False, quoting=csv.QUOTE_ALL)
    print(f"  ✅ Saved {out} with {len(final_df)} rows")
    return out


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2 — Merge duplicate GENE entries (case-insensitive, per FILE+SEQUENCE)
# ─────────────────────────────────────────────────────────────────────────────

def step2_case_merge(input_csv):
    print("\n── Step 2: Case-insensitive gene deduplication ──")
    df = pd.read_csv(input_csv)

    df["_merge_key"] = (
        df["#FILE"].str.strip() + "||" +
        df["SEQUENCE"].str.strip() + "||" +
        df["GENE"].str.strip().str.upper()
    )

    db_merged = (
        df.groupby("_merge_key")["DATABASE"]
        .apply(merge_db_tokens)
        .reset_index()
        .rename(columns={"DATABASE": "DATABASE_merged"})
    )

    df["BEST_DB_RANK"] = df["DATABASE"].map(best_db_rank)
    df_best = (
        df.sort_values("BEST_DB_RANK", ascending=True)
        .drop_duplicates(subset="_merge_key", keep="first")
        .copy()
    )

    df_final = df_best.merge(db_merged, on="_merge_key", how="left")
    df_final["DATABASE"] = df_final["DATABASE_merged"]
    df_final["GENE"] = df_final["GENE"].str.upper()
    df_final = df_final.drop(columns=["_merge_key", "DATABASE_merged"]).reset_index(drop=True)

    print(f"  Before: {len(df)} rows  →  After: {len(df_final)} rows  (removed {len(df) - len(df_final)})")
    out = "deduplicated_ARGs_rare_casemerged.csv"
    df_final.to_csv(out, index=False)
    print(f"  ✅ Saved {out}")
    return out


# ─────────────────────────────────────────────────────────────────────────────
# STEP 3 — Fill missing RESISTANCE via gene-level lookup
# ─────────────────────────────────────────────────────────────────────────────

def step3_fill_resistance(input_csv):
    print("\n── Step 3: Fill missing RESISTANCE annotations ──")
    df = pd.read_csv(input_csv)

    print(f"  With RESISTANCE: {df['RESISTANCE'].notna().sum()}  |  Without: {df['RESISTANCE'].isna().sum()}")

    df["_gene_key"] = df["GENE"].str.strip().str.upper()

    lookup = (
        df.groupby("_gene_key")["RESISTANCE"]
        .apply(merge_resistance_tokens)
        .reset_index()
        .rename(columns={"RESISTANCE": "RESISTANCE_filled"})
    )

    df = df.merge(lookup, on="_gene_key", how="left")
    df["RESISTANCE"] = df["RESISTANCE"].fillna(df["RESISTANCE_filled"])
    df = df.drop(columns=["_gene_key", "RESISTANCE_filled"]).reset_index(drop=True)

    print(f"  With RESISTANCE: {df['RESISTANCE'].notna().sum()}  |  Without: {df['RESISTANCE'].isna().sum()}")
    out = "deduplicated_ARGs_rare_resistance_filled.csv"
    df.to_csv(out, index=False)
    print(f"  ✅ Saved {out}")
    return out


# ─────────────────────────────────────────────────────────────────────────────
# STEP 4 — Manual gene remapping + direct value injection
# ─────────────────────────────────────────────────────────────────────────────

def step4_manual_overrides(input_csv):
    print("\n── Step 4: Manual resistance overrides ──")
    df = pd.read_csv(input_csv)

    print(f"  With RESISTANCE: {df['RESISTANCE'].notna().sum()}  |  Without: {df['RESISTANCE'].isna().sum()}")

    df["_gene_key"] = df["GENE"].str.strip().str.upper()

    # Remap gene keys so they inherit resistance from their source gene
    for target, source in MANUAL_MAPPING.items():
        df.loc[df["_gene_key"] == target.upper(), "_gene_key"] = source.upper()

    lookup = (
        df.groupby("_gene_key")["RESISTANCE"]
        .apply(merge_resistance_tokens)
        .reset_index()
        .rename(columns={"RESISTANCE": "RESISTANCE_filled"})
    )

    # Inject direct values into lookup
    for gene, res_value in MANUAL_VALUES.items():
        mask = lookup["_gene_key"] == gene.upper()
        if mask.any():
            lookup.loc[mask, "RESISTANCE_filled"] = res_value
        else:
            lookup = pd.concat(
                [lookup, pd.DataFrame([{"_gene_key": gene.upper(), "RESISTANCE_filled": res_value}])],
                ignore_index=True,
            )

    df = df.merge(lookup, on="_gene_key", how="left")
    df["RESISTANCE"] = df["RESISTANCE"].replace(["None", "none", "nan"], pd.NA)
    df["RESISTANCE"] = df["RESISTANCE"].fillna(df["RESISTANCE_filled"])
    df = df.drop(columns=["_gene_key", "RESISTANCE_filled"]).reset_index(drop=True)

    print(f"  With RESISTANCE: {df['RESISTANCE'].notna().sum()}  |  Without: {df['RESISTANCE'].isna().sum()}")
    out = "deduplicated_ARGs_rare_resistance_filled_final.csv"
    df.to_csv(out, index=False)
    print(f"  ✅ Saved {out}")
    return out


# ─────────────────────────────────────────────────────────────────────────────
# STEP 5+6 — Map labels → drug classes, classify & explode
# ─────────────────────────────────────────────────────────────────────────────

def step5_map_and_classify(input_csv):
    print("\n── Steps 5+6: Map resistance labels & classify ──")
    df = pd.read_csv(input_csv)

    df["RESISTANCE_MAPPED"] = df["RESISTANCE"].apply(map_resistance_list)

    rows = []
    for _, row in df.iterrows():
        mapped = row["RESISTANCE_MAPPED"]
        classification, classes = classify_resistance(mapped)

        if classification == "multidrug":
            new_row = row.copy()
            new_row["class"] = "multidrug"
            rows.append(new_row)
        elif len(classes) == 2:
            for c in classes:
                new_row = row.copy()
                new_row["class"] = c
                rows.append(new_row)
        elif len(classes) == 1:
            new_row = row.copy()
            new_row["class"] = classes[0]
            rows.append(new_row)
        else:
            new_row = row.copy()
            new_row["class"] = ""
            rows.append(new_row)

    df_final = pd.DataFrame(rows).reset_index(drop=True)
    out = "deduplicated_ARGs_rare_mapped.csv"
    df_final.to_csv(out, index=False)
    print(f"  ✅ Saved {out} with {len(df_final)} rows (from {len(df)} input rows)")
    return out


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    f1 = step1_deduplicate(INPUT_TAB)
    f2 = step2_case_merge(f1)
    f3 = step3_fill_resistance(f2)
    f4 = step4_manual_overrides(f3)
    f5 = step5_map_and_classify(f4)
    print(f"\n🎉 Pipeline complete. Final output: {f5}")
