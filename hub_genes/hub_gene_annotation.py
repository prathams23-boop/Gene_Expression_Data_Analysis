#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python3
"""
annotate_hub_genes.py

Takes a hub-gene CSV (columns: Module, EnsemblID, GeneSymbol, kME) and builds a
consolidated annotation table by querying public APIs:

  - MyGene.info  -> official name, NCBI/RefSeq functional summary, aliases, Entrez ID
  - cBioPortal   -> % of TCGA-BRCA patients with a mutation / copy-number alteration

Survival (KM Plotter) and subtype-expression (UALCAN/GEPIA2) plots are NOT scriptable;
those are web-gated. Use the output CSV to pick the top 2-3 genes per module, then
enter just those into those sites by hand.

Default run needs no arguments:
    python annotate_hub_genes.py
"""

import argparse
import sys
import time
import re
import pandas as pd
import requests
import urllib3

MYGENE_URL = "https://mygene.info/v3/query"
CALLER = "hub_gene_annotation_script"

# TLS verification setting applied to every request.
#   True        -> verify normally (default)
#   False       -> skip verification (set by --insecure; needed behind a
#                  TLS-inspecting campus/corporate proxy presenting its own cert)
#   "<path>"    -> verify against a custom CA bundle (set by --ca-bundle)
VERIFY = True


def strip_version(ensembl_id):
    """ENSG00000141510.14 -> ENSG00000141510"""
    return re.sub(r"\..*$", "", str(ensembl_id))


def load_hub_csv(path):
    df = pd.read_csv(path)
    required = {"Module", "EnsemblID", "GeneSymbol", "kME"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: input CSV is missing columns: {sorted(missing)}")
    df["ens_clean"] = df["EnsemblID"].map(strip_version)
    return df


# ---------------------------------------------------------------------------
# MyGene.info: biological summaries, keyed by Ensembl gene ID (robust vs symbols)
# ---------------------------------------------------------------------------
def fetch_mygene(ens_ids, symbols_by_ens):
    """
    Returns dict keyed by clean Ensembl ID:
        {ens: {"name":..., "summary":..., "alias":..., "entrez":...}}
    Queries by Ensembl ID first; falls back to symbol for any that miss.
    """
    out = {}

    def _post(query_list, scopes):
        payload = {
            "q": ",".join(query_list),
            "scopes": scopes,
            "fields": "symbol,name,summary,alias,entrezgene",
            "species": "human",
        }
        r = requests.post(MYGENE_URL, data=payload, timeout=60, verify=VERIFY)
        r.raise_for_status()
        return r.json()

    # Pass 1: query by Ensembl gene ID
    hits = _post(ens_ids, "ensembl.gene")
    for hit in hits:
        if hit.get("notfound"):
            continue
        q = hit.get("query")  # the ensembl id we sent
        if q in out:
            continue  # keep first (best-scoring) hit per query
        alias = hit.get("alias")
        if isinstance(alias, list):
            alias = "; ".join(alias)
        out[q] = {
            "name": hit.get("name", ""),
            "summary": hit.get("summary", ""),
            "alias": alias or "",
            "entrez": hit.get("entrezgene", ""),
        }

    # Pass 2: fall back to symbol for anything not resolved
    unresolved = [e for e in ens_ids if e not in out]
    if unresolved:
        sym_to_ens = {}
        for e in unresolved:
            s = symbols_by_ens.get(e)
            if s:
                sym_to_ens.setdefault(s, e)
        if sym_to_ens:
            hits2 = _post(list(sym_to_ens.keys()), "symbol")
            for hit in hits2:
                if hit.get("notfound"):
                    continue
                q = hit.get("query")  # the symbol we sent
                e = sym_to_ens.get(q)
                if not e or e in out:
                    continue
                alias = hit.get("alias")
                if isinstance(alias, list):
                    alias = "; ".join(alias)
                out[e] = {
                    "name": hit.get("name", ""),
                    "summary": hit.get("summary", ""),
                    "alias": alias or "",
                    "entrez": hit.get("entrezgene", ""),
                }
    return out


# ---------------------------------------------------------------------------
# cBioPortal: mutation + copy-number alteration frequency in TCGA-BRCA
# ---------------------------------------------------------------------------
def fetch_cbioportal(entrez_ids, base, study):
    """
    Returns dict keyed by entrez id (int):
        {entrez: {"mut_pct":float, "cna_pct":float, "altered_pct":float}}
    Two batch POSTs (mutations, CNA) + two sample-list lookups. Any failure
    raises; caller handles it and leaves these columns blank.
    """
    base = base.rstrip("/")
    mut_profile = f"{study}_mutations"
    cna_profile = f"{study}_gistic"
    seq_list = f"{study}_sequenced"
    cna_list = f"{study}_cna"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    entrez_ids = [int(e) for e in entrez_ids]

    def sample_ids(list_id):
        r = requests.get(f"{base}/sample-lists/{list_id}", timeout=60, verify=VERIFY)
        r.raise_for_status()
        return set(r.json().get("sampleIds", []))

    seq_samples = sample_ids(seq_list)
    cna_samples = sample_ids(cna_list)

    # Mutations
    r = requests.post(
        f"{base}/molecular-profiles/{mut_profile}/mutations/fetch",
        params={"projection": "ID"},
        json={"sampleListId": seq_list, "entrezGeneIds": entrez_ids},
        headers=headers, timeout=120, verify=VERIFY,
    )
    r.raise_for_status()
    mut_samples = {}  # entrez -> set(sampleId)
    for rec in r.json():
        g = rec.get("entrezGeneId")
        s = rec.get("sampleId")
        if g is not None and s is not None:
            mut_samples.setdefault(g, set()).add(s)

    # Copy-number (deep deletions + amplifications only)
    r = requests.post(
        f"{base}/molecular-profiles/{cna_profile}/discrete-copy-number/fetch",
        params={"discreteCopyNumberEventType": "HOMDEL_AND_AMPL", "projection": "ID"},
        json={"sampleListId": cna_list, "entrezGeneIds": entrez_ids},
        headers=headers, timeout=120, verify=VERIFY,
    )
    r.raise_for_status()
    cna_alt = {}  # entrez -> set(sampleId)
    for rec in r.json():
        g = rec.get("entrezGeneId")
        s = rec.get("sampleId")
        if g is not None and s is not None:
            cna_alt.setdefault(g, set()).add(s)

    n_seq = max(len(seq_samples), 1)
    n_cna = max(len(cna_samples), 1)
    profiled_union = seq_samples | cna_samples
    n_union = max(len(profiled_union), 1)

    out = {}
    for g in entrez_ids:
        ms = mut_samples.get(g, set())
        cs = cna_alt.get(g, set())
        out[g] = {
            "mut_pct": round(100.0 * len(ms) / n_seq, 1),
            "cna_pct": round(100.0 * len(cs) / n_cna, 1),
            "altered_pct": round(100.0 * len(ms | cs) / n_union, 1),
        }
    return out


def main():
    ap = argparse.ArgumentParser(description="Consolidated hub-gene annotation CSV.")
    ap.add_argument("--input", default="hub_genes_signed_top20.csv",
                    help="Input hub-gene CSV (Module, EnsemblID, GeneSymbol, kME).")
    ap.add_argument("--output", default="hub_genes_annotated.csv",
                    help="Output consolidated CSV.")
    ap.add_argument("--no-cbioportal", action="store_true",
                    help="Skip the cBioPortal alteration-frequency lookup.")
    ap.add_argument("--cbioportal-base", default="https://www.cbioportal.org/api")
    ap.add_argument("--cbioportal-study", default="brca_tcga_pan_can_atlas_2018")
    ap.add_argument("--insecure", action="store_true",
                    help="Skip TLS certificate verification. Use this behind a "
                         "TLS-inspecting campus/corporate proxy. Safe here because "
                         "these APIs serve only public, read-only gene data.")
    ap.add_argument("--ca-bundle", default=None,
                    help="Path to a custom CA certificate bundle (the proper fix "
                         "if your IT team provides the institutional root cert).")
    args = ap.parse_args()

    global VERIFY
    if args.ca_bundle:
        VERIFY = args.ca_bundle
    elif args.insecure:
        VERIFY = False
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
        print("WARNING: TLS certificate verification is OFF (--insecure).")

    df = load_hub_csv(args.input)
    n_modules = df["Module"].nunique()
    uniq_ens = df["ens_clean"].drop_duplicates().tolist()
    symbols_by_ens = dict(zip(df["ens_clean"], df["GeneSymbol"]))
    print(f"Loaded {len(df)} rows, {len(uniq_ens)} unique genes, {n_modules} modules.")

    # --- MyGene summaries ---
    print("Querying MyGene.info for gene summaries ...")
    mg = fetch_mygene(uniq_ens, symbols_by_ens)
    print(f"  resolved {len(mg)} / {len(uniq_ens)} genes.")

    df["Full_Name"] = df["ens_clean"].map(lambda e: mg.get(e, {}).get("name", ""))
    df["Aliases"] = df["ens_clean"].map(lambda e: mg.get(e, {}).get("alias", ""))
    df["Functional_Summary"] = df["ens_clean"].map(
        lambda e: mg.get(e, {}).get("summary", "") or "No summary available.")
    df["EntrezID"] = df["ens_clean"].map(lambda e: mg.get(e, {}).get("entrez", ""))

    # --- cBioPortal alteration frequencies ---
    df["BRCA_Mutation_pct"] = ""
    df["BRCA_CNA_pct"] = ""
    df["BRCA_Altered_pct"] = ""
    if not args.no_cbioportal:
        entrez_ids = sorted({int(v["entrez"]) for v in mg.values()
                             if str(v.get("entrez", "")).strip().isdigit()})
        if entrez_ids:
            print(f"Querying cBioPortal ({args.cbioportal_study}) for "
                  f"{len(entrez_ids)} genes ...")
            try:
                cb = fetch_cbioportal(entrez_ids, args.cbioportal_base,
                                      args.cbioportal_study)
                ens_to_entrez = {e: mg[e].get("entrez") for e in mg}
                def _get(e, key):
                    ez = ens_to_entrez.get(e)
                    if ez is None or str(ez) == "":
                        return ""
                    return cb.get(int(ez), {}).get(key, "")
                df["BRCA_Mutation_pct"] = df["ens_clean"].map(lambda e: _get(e, "mut_pct"))
                df["BRCA_CNA_pct"] = df["ens_clean"].map(lambda e: _get(e, "cna_pct"))
                df["BRCA_Altered_pct"] = df["ens_clean"].map(lambda e: _get(e, "altered_pct"))
                print("  cBioPortal lookup complete.")
            except Exception as exc:
                print(f"  WARNING: cBioPortal lookup failed ({exc}). "
                      f"Leaving those columns blank.")

    cols = ["Module", "GeneSymbol", "EnsemblID", "kME", "Full_Name", "Aliases",
            "EntrezID", "BRCA_Mutation_pct", "BRCA_CNA_pct", "BRCA_Altered_pct",
            "Functional_Summary"]
    df[cols].to_csv(args.output, index=False)
    print(f"Wrote {args.output} ({len(df)} rows, {len(cols)} columns).")


if __name__ == "__main__":
    main()

