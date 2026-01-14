import os
import time
import pandas as pd
import requests
from tqdm import tqdm

FIND_URL = "https://rest.kegg.jp/find/genes/{query}"
AASEQ_URL = "https://rest.kegg.jp/get/{entry}/aaseq"

INPUT_FILE = "/mnt/NFS/fengch/PAPER/drug/glu_ess.xlsx"  
COLUMN_NAME = "gene"  
OUTPUT_FASTA = "/mnt/NFS/fengch/PAPER/drug/glu_ess.faa"  
LOG_CSV = "/mnt/NFS/fengch/PAPER/drug/kegg_fetch_log.csv"  
SLEEP_TIME = 0.2  

def smart_read_ids(path, colname=None):
    ext = os.path.splitext(path)[1].lower()
    if ext in [".xlsx", ".xls"]:
        df = pd.read_excel(path)
    elif ext == ".csv":
        df = pd.read_csv(path)
    elif ext in [".tsv", ".txt"]:
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception:
            df = pd.read_csv(path)
    else:
        raise ValueError(f"failure: {ext}")
    s = df[colname] if (colname and colname in df.columns) else df.iloc[:,0]
    ids = (s.astype(str).str.strip()).dropna()
    return [x for x in ids if x and x.lower() != "nan"]

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
    txt = http_get(FIND_URL.format(query=cnag))
    if not txt:
        return None, []
    rows = [ln for ln in txt.strip().splitlines() if ln]
    parsed = []
    for ln in rows:
        left, right = ln.split("\t", 1)
        parsed.append((left.strip(), right.strip()))
    exact = [p for p in parsed if p[0].split(":")[-1].upper() == cnag.upper()]
    if exact:
        return exact[0][0], parsed 
    return parsed[0][0], parsed

def kegg_get_aaseq(entry):
    txt = http_get(AASEQ_URL.format(entry=entry))
    if not txt:
        return None, None
    blocks = txt.strip().split("\n>")
    first = blocks[0] if txt.startswith(">") else ">" + blocks[0]
    lines = first.strip().splitlines()
    header = lines[0].lstrip(">").strip()
    seq = "".join(x.strip() for x in lines[1:])
    return header, seq

def main():
    ids = smart_read_ids(INPUT_FILE, COLUMN_NAME)
    print(f"{len(ids)} genes")

    os.makedirs(os.path.dirname(os.path.abspath(OUTPUT_FASTA)) or ".", exist_ok=True)
    logs = []
    ok = 0
    with open(OUTPUT_FASTA, "w", encoding="utf-8") as fw:
        for cnag in tqdm(ids, desc="Fetching KEGG AAseq"):
            try:
                best_entry, candidates = kegg_find_genes(cnag)
                if not best_entry:
                    logs.append({"gene": cnag, "status": "NOT_FOUND_IN_FIND", "best_entry": None, "candidates": ""})
                else:
                    header, seq = kegg_get_aaseq(best_entry)
                    if header and seq:
                        fw.write(f">{header}\n")
                        for i in range(0, len(seq), 60):
                            fw.write(seq[i:i+60] + "\n")
                        ok += 1
                        logs.append({
                            "gene": cnag,
                            "status": "OK",
                            "best_entry": best_entry,
                            "cand_count": len(candidates),
                            "length": len(seq)
                        })
                    else:
                        logs.append({
                            "gene": cnag,
                            "status": "NO_AASEQ",
                            "best_entry": best_entry,
                            "cand_count": len(candidates)
                        })
            except Exception as e:
                logs.append({"gene": cnag, "status": f"ERROR: {e}"})
            time.sleep(SLEEP_TIME)

    pd.DataFrame(logs).to_csv(LOG_CSV, index=False)
    fail = len(ids) - ok
    print(f"done {ok}，failure {fail}。FASTA -> {OUTPUT_FASTA}；log -> {LOG_CSV}")
    
if __name__ == "__main__":
    main()
