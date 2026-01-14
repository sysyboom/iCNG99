import os
import time
import pandas as pd
import requests
from tqdm import tqdm

FIND_URL = "https://rest.kegg.jp/find/genes/{query}"
AASEQ_URL = "https://rest.kegg.jp/get/{entry}/aaseq"

STATUS_FILE   = "/mnt/NFS/fengch/new/models/paper/gene_essentiality_comparison.xlsx"  
SHEET_NAME    = 0          
GENE_COL      = "Gene_ID"  
MODEL_COL     = "model_predicted"
EXPT_COL      = "experimental_essential"

OUTPUT_FASTA  = "/mnt/NFS/fengch/new/models/paper/FN_kegg.faa"       
LOG_CSV       = "/mnt/NFS/fengch/PAPER/FN_kegg_log.csv"   
SLEEP_TIME    = 0.2       
TIMEOUT       = 30
RETRIES       = 2
BACKOFF       = 1.6

def smart_read_table(path, sheet=None):
    ext = os.path.splitext(path)[1].lower()
    if ext in [".xlsx", ".xls"]:
        return pd.read_excel(path, sheet_name=sheet)
    elif ext == ".csv":
        return pd.read_csv(path)
    elif ext in [".tsv", ".txt"]: 
        try:
            return pd.read_csv(path, sep="\t")
        except Exception:
            return pd.read_csv(path)
    else:
        raise ValueError(f"false: {ext}")


def http_get(url, timeout=30, retries=2, backoff=1.6):
    last_err = None
    for i in range(retries + 1):
        try:
            r = requests.get(url, timeout=timeout)
            if r.status_code == 200 and r.text.strip():
                return r.text
            if r.status_code == 404:
                return None
            last_err = f"HTTP {r.status_code}"
        except Exception as e:
            last_err = str(e)
        time.sleep(backoff ** i)
    if "HTTP 404" in str(last_err):
        return None
    raise RuntimeError(f"failure: {url} | {last_err}")


def kegg_find_genes(cnag):
    txt = http_get(FIND_URL.format(query=cnag), timeout=TIMEOUT, retries=RETRIES, backoff=BACKOFF)
    if not txt:
        return None, []
    rows = [ln for ln in txt.strip().splitlines() if ln]
    parsed = []
    for ln in rows:
        if "\t" in ln:
            left, right = ln.split("\t", 1)
            parsed.append((left.strip(), right.strip()))
    exact = [p for p in parsed if p[0].split(":")[-1].upper() == cnag.upper()]
    if exact:
        return exact[0][0], parsed
    return parsed[0][0], parsed


def kegg_get_aaseq(entry):
    txt = http_get(AASEQ_URL.format(entry=entry), timeout=TIMEOUT, retries=RETRIES, backoff=BACKOFF)
    if not txt:
        return None, None
    blocks = txt.strip().split("\n>")
    first = blocks[0] if txt.startswith(">") else ">" + blocks[0]
    lines = first.strip().splitlines()
    header = lines[0].lstrip(">").strip()
    seq = "".join(x.strip() for x in lines[1:])
    return header, seq


def write_fasta_record(fh, header, seq, width=60):
    fh.write(f">{header}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i+width] + "\n")


def main():
    df = smart_read_table(STATUS_FILE, SHEET_NAME)
    need_cols = {GENE_COL, MODEL_COL, EXPT_COL}
    if not need_cols.issubset(df.columns):
        raise ValueError(f"lack：{need_cols - set(df.columns)}")

    df[GENE_COL] = df[GENE_COL].astype(str).str.strip()
    df = df.dropna(subset=[GENE_COL])

    fn_genes = (
        df[(df[MODEL_COL] == 0) & (df[EXPT_COL] == 1)][GENE_COL]
        .astype(str).str.strip().str.upper().dropna().unique().tolist()
    )
    print(f"FN: {len(fn_genes)} ")

    os.makedirs(os.path.dirname(os.path.abspath(OUTPUT_FASTA)) or ".", exist_ok=True)
    logs = []
    ok = 0

    with open(OUTPUT_FASTA, "w", encoding="utf-8") as fw:
        for gene in tqdm(fn_genes, desc="Fetching KEGG AAseq (FN)"):
            try:
                best_entry, cands = kegg_find_genes(gene)
                if not best_entry:
                    logs.append({
                        "gene": gene, "status": "NOT_FOUND_IN_FIND",
                        "best_entry": None, "cand_count": 0
                    })
                else:
                    header, seq = kegg_get_aaseq(best_entry)
                    if header and seq:
                        write_fasta_record(fw, header, seq)
                        ok += 1
                        logs.append({
                            "gene": gene, "status": "OK",
                            "best_entry": best_entry,
                            "cand_count": len(cands),
                            "length": len(seq)
                        })
                    else:
                        logs.append({
                            "gene": gene, "status": "NO_AASEQ",
                            "best_entry": best_entry,
                            "cand_count": len(cands)
                        })
            except Exception as e:
                logs.append({"gene": gene, "status": f"ERROR: {e}"})
            time.sleep(SLEEP_TIME)

    pd.DataFrame(logs).to_csv(LOG_CSV, index=False)


if __name__ == "__main__":
    main()
