#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从 CNAG 基因名列表直接在 KEGG 全库查找并提取 AA 序列
- 输入：xlsx/csv/tsv/txt（默认取第一列或用 -c 指定列）
- 逻辑：/find/genes/{CNAG} -> 命中条目(如 cng:CNAG_01234) -> /get/{entry}/aaseq
- 输出：kegg_cnag.faa（蛋白FASTA） + kegg_fetch_log.csv（日志）
依赖：requests, pandas, tqdm   pip install requests pandas tqdm
"""

import os
import time
import pandas as pd
import requests
from tqdm import tqdm

FIND_URL = "https://rest.kegg.jp/find/genes/{query}"
AASEQ_URL = "https://rest.kegg.jp/get/{entry}/aaseq"

# 直接在代码中设置参数
INPUT_FILE = "/mnt/NFS/fengch/PAPER/drug/glu_ess.xlsx"  # 输入文件路径（Excel 文件）
COLUMN_NAME = "gene"  # Excel 中基因名的列名，如果不指定会默认使用第一列
OUTPUT_FASTA = "/mnt/NFS/fengch/PAPER/drug/glu_ess.faa"  # 输出的 FASTA 文件路径
LOG_CSV = "/mnt/NFS/fengch/PAPER/drug/kegg_fetch_log.csv"  # 输出的日志 CSV 文件路径
SLEEP_TIME = 0.2  # 每次请求后的休眠时间，单位为秒

def smart_read_ids(path, colname=None):
    ext = os.path.splitext(path)[1].lower()
    if ext in [".xlsx", ".xls"]:
        df = pd.read_excel(path)
    elif ext == ".csv":
        df = pd.read_csv(path)
    elif ext in [".tsv", ".txt"]:
        # 先尝试tab，再退回逗号
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception:
            df = pd.read_csv(path)
    else:
        raise ValueError(f"不支持的文件格式: {ext}")
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
    raise RuntimeError(f"请求失败: {url} | {last_err}")

def kegg_find_genes(cnag):
    """在 KEGG genes 全库里查找该 CNAG，返回最佳 entry（如 cng:CNAG_01234）"""
    txt = http_get(FIND_URL.format(query=cnag))
    if not txt:
        return None, []
    rows = [ln for ln in txt.strip().splitlines() if ln]
    # 解析 entry 与描述
    parsed = []
    for ln in rows:
        left, right = ln.split("\t", 1)
        parsed.append((left.strip(), right.strip()))
    # 优先 exact：左侧应为 org:CNAG_XXXX 且恰好匹配
    exact = [p for p in parsed if p[0].split(":")[-1].upper() == cnag.upper()]
    if exact:
        return exact[0][0], parsed  # 返回最佳和所有候选
    # 否则取第一条（通常是最相关）
    return parsed[0][0], parsed

def kegg_get_aaseq(entry):
    txt = http_get(AASEQ_URL.format(entry=entry))
    if not txt:
        return None, None
    # 可能多条，这里取第一条
    blocks = txt.strip().split("\n>")
    first = blocks[0] if txt.startswith(">") else ">" + blocks[0]
    lines = first.strip().splitlines()
    header = lines[0].lstrip(">").strip()
    seq = "".join(x.strip() for x in lines[1:])
    return header, seq

def main():
    ids = smart_read_ids(INPUT_FILE, COLUMN_NAME)
    print(f"读取到 {len(ids)} 个 CNAG 基因")

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
    print(f"完成：成功 {ok}，失败 {fail}。FASTA -> {OUTPUT_FASTA}；日志 -> {LOG_CSV}")

if __name__ == "__main__":
    main()
