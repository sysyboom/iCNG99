#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从“基因状态表”中筛出 FN 基因（model_predicted==0 & experimental_essential==1），
然后到 KEGG 拉取对应的蛋白 AA 序列。

输入：
- 一个 Excel/CSV/TSV/TXT 的“基因状态表”，至少包含列：
  Gene_ID, model_predicted, experimental_essential
  （Gene_ID 例如 CNAG_01234）

流程：
- 从状态表筛出 FN 基因 -> 对每个 CNAG 用 KEGG REST:
    /find/genes/{CNAG}  找 entry（如 cng:CNAG_01234）
    /get/{entry}/aaseq  拉 AA 序列
- 输出：
  1) FASTA：FN_kegg.faa
  2) 日志：FN_kegg_log.csv（gene、entry、长度、失败原因等）
依赖：requests, pandas, tqdm   （pip install requests pandas tqdm）
"""

import os
import time
import pandas as pd
import requests
from tqdm import tqdm

# ---------- KEGG REST ----------
FIND_URL = "https://rest.kegg.jp/find/genes/{query}"
AASEQ_URL = "https://rest.kegg.jp/get/{entry}/aaseq"

# ---------- ☆☆ 需要你改的参数 ☆☆ ----------
STATUS_FILE   = "/mnt/NFS/fengch/new/models/paper/gene_essentiality_comparison.xlsx"  # 含 Gene_ID / model_predicted / experimental_essential
SHEET_NAME    = 0          # Excel 工作表名或索引；CSV/TSV/TXT 无需改
GENE_COL      = "Gene_ID"  # 列名：基因
MODEL_COL     = "model_predicted"
EXPT_COL      = "experimental_essential"

OUTPUT_FASTA  = "/mnt/NFS/fengch/new/models/paper/FN_kegg.faa"       # 输出 FASTA
LOG_CSV       = "/mnt/NFS/fengch/PAPER/FN_kegg_log.csv"   # 日志 CSV
SLEEP_TIME    = 0.2       # 每次请求后的休眠（秒）
TIMEOUT       = 30
RETRIES       = 2
BACKOFF       = 1.6
# -------------------------------------------


def smart_read_table(path, sheet=None):
    ext = os.path.splitext(path)[1].lower()
    if ext in [".xlsx", ".xls"]:
        return pd.read_excel(path, sheet_name=sheet)
    elif ext == ".csv":
        return pd.read_csv(path)
    elif ext in [".tsv", ".txt"]:
        # 优先尝试制表符
        try:
            return pd.read_csv(path, sep="\t")
        except Exception:
            return pd.read_csv(path)
    else:
        raise ValueError(f"不支持的文件格式: {ext}")


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
    raise RuntimeError(f"请求失败: {url} | {last_err}")


def kegg_find_genes(cnag):
    """
    在 KEGG genes 里查 CNAG，返回 (最佳entry, 候选列表)。
    最佳优先 exact 命中（左侧 org:CNAG_xxx 与查询一致），否则取第一条。
    """
    txt = http_get(FIND_URL.format(query=cnag), timeout=TIMEOUT, retries=RETRIES, backoff=BACKOFF)
    if not txt:
        return None, []
    rows = [ln for ln in txt.strip().splitlines() if ln]
    parsed = []
    for ln in rows:
        if "\t" in ln:
            left, right = ln.split("\t", 1)
            parsed.append((left.strip(), right.strip()))
    # exact 命中
    exact = [p for p in parsed if p[0].split(":")[-1].upper() == cnag.upper()]
    if exact:
        return exact[0][0], parsed
    return parsed[0][0], parsed


def kegg_get_aaseq(entry):
    """
    拉取 AA 序列，返回 (header, seq)。如果返回多条，这里取第一条。
    """
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
    # 读取状态表并筛出 FN
    df = smart_read_table(STATUS_FILE, SHEET_NAME)
    need_cols = {GENE_COL, MODEL_COL, EXPT_COL}
    if not need_cols.issubset(df.columns):
        raise ValueError(f"输入文件缺少列：{need_cols - set(df.columns)}")

    # 去空白、标准化
    df[GENE_COL] = df[GENE_COL].astype(str).str.strip()
    df = df.dropna(subset=[GENE_COL])

    # FN：模型0 & 实验1
    fn_genes = (
        df[(df[MODEL_COL] == 0) & (df[EXPT_COL] == 1)][GENE_COL]
        .astype(str).str.strip().str.upper().dropna().unique().tolist()
    )
    print(f"[INFO] 筛到 FN 基因 {len(fn_genes)} 个。")

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
    print(f"[DONE] 成功 {ok} / {len(fn_genes)}。FASTA -> {OUTPUT_FASTA}；日志 -> {LOG_CSV}")


if __name__ == "__main__":
    main()
