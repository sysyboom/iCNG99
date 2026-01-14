import os
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

def export_up_down_logFC_with_padj(
    in_file: str,
    out_up_xlsx: str,
    out_down_xlsx: str,
    lfc_threshold: float = 0.858,
    padj_threshold: float = 0.05,
    pseudocount: float = 1e-3,
    sheet_name = 0
):

    ext = os.path.splitext(in_file)[1].lower()
    if ext in [".xlsx", ".xls"]:
        df = pd.read_excel(in_file, sheet_name=sheet_name)
    else:
        sep = "\t" if ext in [".tsv", ".txt"] else ","
        df = pd.read_csv(in_file, sep=sep)

    if df.shape[1] < 5:
        raise ValueError("Need at least 1 ID column + >=4 sample columns (>=2 control and >=2 treatment).")

    reaction = df.iloc[:, 0].astype(str)

    X = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    n_samples = X.shape[1]
    if n_samples % 2 != 0:
        raise ValueError(f"Sample columns must be even (split half/half). Got {n_samples} sample columns.")

    mid = n_samples // 2
    control = X[:, :mid]
    treat = X[:, mid:]

    c_mean = np.nanmean(control, axis=1)
    t_mean = np.nanmean(treat, axis=1)
    log2fc = np.log2((t_mean + pseudocount) / (c_mean + pseudocount))

    control_log = np.log2(control + pseudocount)
    treat_log = np.log2(treat + pseudocount)

    pvals = np.array([
        ttest_ind(treat_log[i, :], control_log[i, :], equal_var=False, nan_policy="omit").pvalue
        for i in range(X.shape[0])
    ], dtype=float)

    ok = np.isfinite(pvals)
    padj = np.full_like(pvals, np.nan, dtype=float)

    p_ok = pvals[ok]
    m = p_ok.size
    if m > 0:
        order = np.argsort(p_ok)
        ranked_p = p_ok[order]
        q = ranked_p * m / (np.arange(1, m + 1))
        q = np.minimum.accumulate(q[::-1])[::-1]
        q = np.clip(q, 0, 1)

        padj_ok = np.empty_like(q)
        padj_ok[order] = q
        padj[ok] = padj_ok

    res = pd.DataFrame({
        "reaction": reaction,
        "value": log2fc,
        "padj": padj
    }).dropna(subset=["reaction", "value", "padj"])

    res = res[res["padj"] <= float(padj_threshold)]

    up = res[res["value"] >= float(lfc_threshold)].sort_values("value", ascending=False)
    down = res[res["value"] <= -float(lfc_threshold)].sort_values("value", ascending=True)

    up.to_excel(out_up_xlsx, index=False)
    down.to_excel(out_down_xlsx, index=False)

    print(f"Filtered Up-regulated reactions: {up.shape[0]}")
    print(f"Filtered Down-regulated reactions: {down.shape[0]}")

    return up.shape[0], down.shape[0]


if __name__ == "__main__":
    export_up_down_logFC_with_padj(
        in_file="/mnt/NFS/fengch/TPM/vivo/vivo_reactions.xlsx",    
        out_up_xlsx="/mnt/NFS/fengch/TPM/vivo/vivo585_up.xlsx",
        out_down_xlsx="/mnt/NFS/fengch/TPM/vivo/vivo585_down.xlsx",
        lfc_threshold=0.585,
        padj_threshold=0.05,
        pseudocount=1e-3,
        sheet_name=0
    )
