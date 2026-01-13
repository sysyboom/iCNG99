#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import os
import pandas as pd
from cobra.io import read_sbml_model
from cobra.flux_analysis import flux_variability_analysis

# ======================
# 参数（按需修改）
# ======================
MODEL_CONTROL = "/mnt/NFS/fengch/new/models/paper/merge_after_YPD.xml"
MODEL_TREATED = "/mnt/NFS/fengch/new/models/paper/merge_after_vivo.xml"

UP_XLSX    = "/mnt/NFS/fengch/TPM/vivo/vivo585_up.xlsx"
DOWN_XLSX  = "/mnt/NFS/fengch/TPM/vivo/vivo585_down.xlsx"

OUT_UP_IDX   = "/mnt/NFS/fengch/TPM/vivo/vivo585_up_index.xlsx"
OUT_UP_VAL   = "/mnt/NFS/fengch/TPM/vivo/vivo585_up_value.xlsx"
OUT_DN_IDX   = "/mnt/NFS/fengch/TPM/vivo/vivo585_down_index.xlsx"
OUT_DN_VAL   = "/mnt/NFS/fengch/TPM/vivo/vivo585_down_value.xlsx"

FVA_FRACTION = 0.99     # FVA 的 fraction_of_optimum（0.95 可更快）
TOL_ZERO     = 1e-6     # 判定“0通量”容差
Q_DMID       = 0.7     # Δmid 分位数阈值（0.80 更松，0.90 更紧）
Q_DWIDTH     = 0.7     # Δwidth 分位数阈值
OVERLAP_MAX  = 0.6     # 区间重叠比上限（0.40 更严格）
TOL_NONZERO  = 1e-6     # 非零保护阈值：至少一侧远离 0

# 仅对 up/down 涉及的反应做 FVA（显著提速）
FVA_ON_SUBSET = True
FVA_PROCESSES = 8       # 并行进程数（按你的CPU改）

# 可选：排除 EX_/DM_/SK_ 等交换/需求反应（非研究重点时建议 True）
EXCLUDE_EX_DM_SK = False
# ======================

def _normalize_fva_cols(df: pd.DataFrame) -> pd.DataFrame:
    """把 FVA 结果列名统一为 minimum/maximum（兼容 min/max）"""
    rename_map = {}
    for c in df.columns:
        lc = c.lower()
        if lc in ("min", "minimum"):
            rename_map[c] = "minimum"
        elif lc in ("max", "maximum"):
            rename_map[c] = "maximum"
    if rename_map:
        df = df.rename(columns=rename_map)
    return df

def main():
    # ---------- 1. 读模型 ----------
    model_control = read_sbml_model(MODEL_CONTROL)
    model_treated = read_sbml_model(MODEL_TREATED)

    # ---------- 2. 读上下调表 ----------
    def safe_read_excel(path):
        if path is None or not os.path.exists(path):
            print(f"⚠️ 未找到文件: {path} → 视为空表")
            return pd.DataFrame(columns=["reaction", "value"])
        df = pd.read_excel(path)
        if df.empty:
            print(f"⚠️ 文件为空: {path}")
            return pd.DataFrame(columns=["reaction", "value"])
        return df

    up_df = safe_read_excel(UP_XLSX)
    down_df = safe_read_excel(DOWN_XLSX)

    # ---------- 3. 列检查 ----------
    for name, df in {"up_df": up_df, "down_df": down_df}.items():
        if not all(c in df.columns for c in ["reaction", "value"]):
            print(f"⚠️ {name} 缺列 → 自动补充空列")
            for c in ["reaction", "value"]:
                if c not in df.columns:
                    df[c] = []

    # 可选：排除 EX_/DM_/SK_ 类反应
    if EXCLUDE_EX_DM_SK:
        pat = re.compile(r"^(EX_|DM_|SK_)")
        up_df   = up_df  [~up_df["reaction"].astype(str).str.match(pat)].copy()
        down_df = down_df[~down_df["reaction"].astype(str).str.match(pat)].copy()

    print(f"[起始] up={len(up_df)}, down={len(down_df)}")

    # ---------- 4. 两模型 FVA ----------
    rxn_list = None
    if FVA_ON_SUBSET:
        # 仅在 up/down 出现过的反应上做 FVA
        model_rxn_ids = {r.id for r in model_control.reactions}
        rxn_list = [r for r in set(up_df["reaction"]).union(set(down_df["reaction"])) if r in model_rxn_ids]
        print(f"[FVA范围] 仅在 up/down 涉及的 {len(rxn_list)} 条反应上做 FVA")

    fva_control = flux_variability_analysis(
        model_control,
        reaction_list=rxn_list,
        fraction_of_optimum=FVA_FRACTION,
        processes=FVA_PROCESSES
    )
    fva_treated = flux_variability_analysis(
        model_treated,
        reaction_list=rxn_list,
        fraction_of_optimum=FVA_FRACTION,
        processes=FVA_PROCESSES
    )

    fva_control = _normalize_fva_cols(fva_control)
    fva_treated = _normalize_fva_cols(fva_treated)

    # ---------- 5. 两条件都完全不活跃 ----------
    inactive_control = set(fva_control[(fva_control["minimum"].abs() < TOL_ZERO) &
                                       (fva_control["maximum"].abs() < TOL_ZERO)].index)
    inactive_treated = set(fva_treated[(fva_treated["minimum"].abs() < TOL_ZERO) &
                                       (fva_treated["maximum"].abs() < TOL_ZERO)].index)
    inactive_both = inactive_control & inactive_treated

    # ---------- 6. 去掉完全不活跃 ----------
    filtered_up_df   = up_df  [~up_df["reaction"].isin(inactive_both)].copy()
    filtered_down_df = down_df[~down_df["reaction"].isin(inactive_both)].copy()
    print(f"[剔除两条件皆不活跃] up={len(filtered_up_df)}, down={len(filtered_down_df)}")

    # ---------- 7. 仅用 FVA 做主筛（Δmid/Δwidth + 重叠比 + 非零保护） ----------
    common = fva_control.index.intersection(fva_treated.index)
    loA, hiA = fva_control.loc[common, "minimum"], fva_control.loc[common, "maximum"]
    loB, hiB = fva_treated.loc[common, "minimum"], fva_treated.loc[common, "maximum"]

    # Δmid 与 Δwidth
    dmid   = ((loA + hiA) / 2.0 - (loB + hiB) / 2.0).abs()
    widthA = (hiA - loA).abs()
    widthB = (hiB - loB).abs()
    dwidth = (widthA - widthB).abs()

    cut_dmid   = dmid.quantile(Q_DMID)
    cut_dwidth = dwidth.quantile(Q_DWIDTH)

    # 区间重叠比
    overlap = (pd.concat([hiA, hiB], axis=1).min(axis=1) -
               pd.concat([loA, loB], axis=1).max(axis=1)).clip(lower=0)
    union = (pd.concat([hiA, hiB], axis=1).max(axis=1) -
             pd.concat([loA, loB], axis=1).min(axis=1)).clip(lower=1e-12)
    keep_overlap = (overlap / union) < OVERLAP_MAX

    # 非零保护：至少一侧不全在 0 附近
    farA = pd.concat([loA.abs(), hiA.abs()], axis=1).max(axis=1)
    farB = pd.concat([loB.abs(), hiB.abs()], axis=1).max(axis=1)
    nonzero_either = (farA >= TOL_NONZERO) | (farB >= TOL_NONZERO)

    # 最终保留
    primary = (dmid >= cut_dmid) | (dwidth >= cut_dwidth)
    keep = primary & keep_overlap & nonzero_either
    keep_rxns = set(keep[keep].index)

    filtered_up_df   = filtered_up_df  [filtered_up_df["reaction"].isin(keep_rxns)].copy()
    filtered_down_df = filtered_down_df[filtered_down_df["reaction"].isin(keep_rxns)].copy()

    print(f"[ΔFVA筛选] Δmid≥Q{int(Q_DMID*100)} 或 Δwidth≥Q{int(Q_DWIDTH*100)}，且 overlap<{OVERLAP_MAX}, 非零保护")
    print(f"[结果] up={len(filtered_up_df)}, down={len(filtered_down_df)}")

    # （可选）再收口：Top-K（按 Δmid+Δwidth）
    # 如需启用，解开注释：
    # ENABLE_TOPK = False
    # TOPK = 120
    # if ENABLE_TOPK:
    #     score = (dmid.reindex(common).fillna(0) + dwidth.reindex(common).fillna(0))
    #     union_rxns = pd.Index(filtered_up_df["reaction"]).union(filtered_down_df["reaction"])
    #     topk_rxns = score.reindex(union_rxns).sort_values(ascending=False).head(TOPK).index
    #     filtered_up_df   = filtered_up_df  [filtered_up_df["reaction"].isin(topk_rxns)].copy()
    #     filtered_down_df = filtered_down_df[filtered_down_df["reaction"].isin(topk_rxns)].copy()
    #     print(f"[Top-K] 取前 {TOPK} 条（按 Δmid+Δwidth） → up={len(filtered_up_df)}, down={len(filtered_down_df)}")

    # ---------- 8. 映射反应 ID → 索引（基于 control 模型顺序） ----------
    rxn_to_index = {rxn.id: idx for idx, rxn in enumerate(model_control.reactions)}
    filtered_up_df["index"]   = filtered_up_df["reaction"].map(rxn_to_index)
    filtered_down_df["index"] = filtered_down_df["reaction"].map(rxn_to_index)

    # 丢掉映射失败的行，防止 astype 报错
    filtered_up_df   = filtered_up_df.dropna(subset=["index"])
    filtered_down_df = filtered_down_df.dropna(subset=["index"])

    filtered_up_df["index"]   = filtered_up_df["index"].astype(int) + 1
    filtered_down_df["index"] = filtered_down_df["index"].astype(int) + 1

    # ---------- 9. 输出 ----------
    filtered_up_df[["index"]].to_excel(OUT_UP_IDX, index=False, header=False)
    filtered_up_df[["value"]].to_excel(OUT_UP_VAL, index=False, header=False)

    filtered_down_df[["index"]].to_excel(OUT_DN_IDX, index=False, header=False)
    filtered_down_df[["value"]].to_excel(OUT_DN_VAL, index=False, header=False)

    print("[完成] 已输出四个文件：")
    print(" -", OUT_UP_IDX)
    print(" -", OUT_UP_VAL)
    print(" -", OUT_DN_IDX)
    print(" -", OUT_DN_VAL)

if __name__ == "__main__":
    main()
