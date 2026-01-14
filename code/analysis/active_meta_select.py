import os
import pandas as pd
from cobra.io import read_sbml_model
from cobra.flux_analysis import flux_variability_analysis

MODEL_PATH = "/mnt/NFS/fengch/new/models/paper/merge_after_YPD_heat.xml"
OUT_DIR = "/mnt/NFS/fengch/new_data/new_results_YPD"
os.makedirs(OUT_DIR, exist_ok=True)

FVA_FRACTION = 0.99
FVA_PROCESSES = 8
TOL_ZERO = 1e-9

def main():
    model = read_sbml_model(MODEL_PATH)
    print(f"{len(model.reactions)} reactions, {len(model.metabolites)} metabolites")

    fva_res = flux_variability_analysis(
        model,
        fraction_of_optimum=FVA_FRACTION,
        processes=FVA_PROCESSES
    )

    fva_res = fva_res.rename(columns=lambda c: "minimum" if c.lower().startswith("min") else
                             ("maximum" if c.lower().startswith("max") else c))

    dead_rxns = fva_res[
        (fva_res["minimum"].abs() < TOL_ZERO) &
        (fva_res["maximum"].abs() < TOL_ZERO)
    ].index.tolist()

    active_rxns = [r.id for r in model.reactions if r.id not in dead_rxns]

    active_mets = set()
    for rxn_id in active_rxns:
        rxn = model.reactions.get_by_id(rxn_id)
        for m in rxn.metabolites:
            active_mets.add(m.id)

    active_mets = sorted(list(active_mets))
    df_active = pd.DataFrame({"Active_Metabolite_in_YPD": active_mets})
    df_active.to_csv(os.path.join(OUT_DIR, "YPD_active_metabolites.tsv"),
                     sep="\t", index=False)

    print(f"active metabolites: {len(active_mets)}")

if __name__ == "__main__":
    main()
