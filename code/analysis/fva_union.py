#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from cobra.io import read_sbml_model
from cobra.flux_analysis import flux_variability_analysis

# ======================
#  参数区（按需修改）
# ======================
MODEL_A_PATH = "/mnt/NFS/fengch/new/models/paper/merge_after_YPD_heat.xml"   # 模型1
MODEL_B_PATH = "/mnt/NFS/fengch/new/models/paper/merge_after_YPD_heat.xml"  # 模型2

OUT_DIR = "/mnt/NFS/fengch/TPM/heat"

FVA_FRACTION = 0.8
FVA_PROCESSES = 8
TOL_ZERO = 1e-9         # 判定 active 的阈值

# 输出文件
FVA_A_CSV      = os.path.join(OUT_DIR, "fva_heat.csv")
FVA_B_CSV      = os.path.join(OUT_DIR, "fva_heat.csv")
ACTIVE_ANY_TXT = os.path.join(OUT_DIR, "active_rxns_union_heat.txt")


def run_fva_for_model(model_path, out_csv):
    print(f"读取模型: {model_path}")
    model = read_sbml_model(model_path)
    print(f"模型包含 {len(model.reactions)} 条反应, {len(model.metabolites)} 个代谢物")

    print("开始运行 FVA ...")
    fva_res = flux_variability_analysis(
        model,
        fraction_of_optimum=FVA_FRACTION,
        processes=FVA_PROCESSES
    )  # index 是 reaction.id, 列是 minimum / maximum

    # meta 信息（同样以 ID 为 index）
    rxn_ids, rxn_names, lb_list, ub_list = [], [], [], []
    for rxn in model.reactions:
        rxn_ids.append(rxn.id)
        rxn_names.append(rxn.name)
        lb_list.append(rxn.lower_bound)
        ub_list.append(rxn.upper_bound)

    meta_df = pd.DataFrame({
        "ID": rxn_ids,
        "Name": rxn_names,
        "LB": lb_list,
        "UB": ub_list
    }).set_index("ID")

    # 合并 meta + FVA
    fva_df = meta_df.join(fva_res, how="left")  # index = ID

    # 判定 active
    fva_df["is_active"] = (
        (fva_df["minimum"].abs() > TOL_ZERO) |
        (fva_df["maximum"].abs() > TOL_ZERO)
    )

    num_active = int(fva_df["is_active"].sum())
    print(f"{os.path.basename(model_path)} 中 active 反应数量: {num_active}")

    # 把 index(ID) 变成一列，方便后面用
    fva_df = fva_df.reset_index()  # 现在有一列叫 "ID"

    # 保存 CSV
    fva_df.to_csv(out_csv, index=False)
    print(f"FVA 结果已保存到: {out_csv}")

    return fva_df


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # 模型A
    fva_a = run_fva_for_model(MODEL_A_PATH, FVA_A_CSV)
    fva_a = fva_a.set_index("ID")  # 现在有 ID 列了，这行就不会报错

    # 模型B
    fva_b = run_fva_for_model(MODEL_B_PATH, FVA_B_CSV)
    fva_b = fva_b.set_index("ID")

    # 合并两个模型的 is_active
    all_ids = sorted(set(fva_a.index) | set(fva_b.index))
    union_df = pd.DataFrame(index=all_ids)

    union_df["is_active_A"] = fva_a.reindex(all_ids)["is_active"].fillna(False)
    union_df["is_active_B"] = fva_b.reindex(all_ids)["is_active"].fillna(False)

    union_df["is_active_any"] = union_df["is_active_A"] | union_df["is_active_B"]

    num_any = int(union_df["is_active_any"].sum())
    print(f"两个模型并集中 active 反应数量: {num_any}")

    # 导出“至少在一个模型里活跃”的反应 ID
    active_any_ids = union_df.index[union_df["is_active_any"]].tolist()
    with open(ACTIVE_ANY_TXT, "w") as f:
        for rid in active_any_ids:
            f.write(rid + "\n")
    print(f"active 反应并集列表已保存到: {ACTIVE_ANY_TXT}")


if __name__ == "__main__":
    main()
