#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from cobra.io import read_sbml_model
from cobra.flux_analysis import flux_variability_analysis

MODEL_A_PATH = "/mnt/NFS/fengch/new/models/paper/merge_after_YPD_heat.xml"   
MODEL_B_PATH = "/mnt/NFS/fengch/new/models/paper/merge_after_YPD_heat.xml"  

OUT_DIR = "/mnt/NFS/fengch/TPM/heat"

FVA_FRACTION = 0.8
FVA_PROCESSES = 8
TOL_ZERO = 1e-9        

FVA_A_CSV      = os.path.join(OUT_DIR, "fva_heat.csv")
FVA_B_CSV      = os.path.join(OUT_DIR, "fva_heat.csv")
ACTIVE_ANY_TXT = os.path.join(OUT_DIR, "active_rxns_union_heat.txt")


def run_fva_for_model(model_path, out_csv):
    model = read_sbml_model(model_path)

    fva_res = flux_variability_analysis(
        model,
        fraction_of_optimum=FVA_FRACTION,
        processes=FVA_PROCESSES
    ) 

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

    fva_df = meta_df.join(fva_res, how="left")  # index = ID

    fva_df["is_active"] = (
        (fva_df["minimum"].abs() > TOL_ZERO) |
        (fva_df["maximum"].abs() > TOL_ZERO)
    )

    num_active = int(fva_df["is_active"].sum())
    print(f"{os.path.basename(model_path)} 中 active reactions: {num_active}")

    fva_df = fva_df.reset_index()  
    fva_df.to_csv(out_csv, index=False)
    return fva_df


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    fva_a = run_fva_for_model(MODEL_A_PATH, FVA_A_CSV)
    fva_a = fva_a.set_index("ID")  

    fva_b = run_fva_for_model(MODEL_B_PATH, FVA_B_CSV)
    fva_b = fva_b.set_index("ID")

    all_ids = sorted(set(fva_a.index) | set(fva_b.index))
    union_df = pd.DataFrame(index=all_ids)

    union_df["is_active_A"] = fva_a.reindex(all_ids)["is_active"].fillna(False)
    union_df["is_active_B"] = fva_b.reindex(all_ids)["is_active"].fillna(False)

    union_df["is_active_any"] = union_df["is_active_A"] | union_df["is_active_B"]

    num_any = int(union_df["is_active_any"].sum())
    print(f"all active reactions: {num_any}")

    active_any_ids = union_df.index[union_df["is_active_any"]].tolist()
    with open(ACTIVE_ANY_TXT, "w") as f:
        for rid in active_any_ids:
            f.write(rid + "\n")
    print(f"active reactions done: {ACTIVE_ANY_TXT}")


if __name__ == "__main__":
    main()
