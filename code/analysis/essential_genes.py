import cobra
import pandas as pd

model = cobra.io.read_sbml_model("/mnt/NFS/fengch/TPM/drug/Gal_GIMME.xml")
# reactions_to_modify = {
    #"EX_003": {"lower_bound": -1000, "upper_bound": 1000}, 
    #"EX_111": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_112": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_113": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_114": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_116": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_117": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_118": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_119": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_120": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_121": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_122": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_123": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_124": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_125": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_126": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_127": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_128": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_129": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_130": {"lower_bound": -1000, "upper_bound": 1000},
    #"EX_162": {"lower_bound": -1000, "upper_bound": 1000},
    #"biomass1": {"lower_bound": 5, "upper_bound": 1000},
#}

# for reaction_id, bounds in reactions_to_modify.items():
    # reaction = model.reactions.get_by_id(reaction_id)
    # reaction.lower_bound = bounds["lower_bound"]
    # reaction.upper_bound = bounds["upper_bound"]

new_objective_reaction_id = "biomass1" 
model.objective = new_objective_reaction_id

gene_knockout_results = pd.DataFrame(columns=["Gene", "Wild-type Objective Value", "Knockout Objective Value"])

for gene in model.genes:
    gene_id = gene.id

    knockout_model = model.copy()
    knockout_model.genes.get_by_id(gene_id).knock_out()

    wildtype_solution = model.slim_optimize()
    knockout_solution = knockout_model.slim_optimize()

    gene_knockout_results = gene_knockout_results.append({"Gene": gene_id,
                                                          "Wild-type Objective Value": wildtype_solution,
                                                          "Knockout Objective Value": knockout_solution},
                                                         ignore_index=True)

gene_knockout_results["Objective Value Change"] = gene_knockout_results["Wild-type Objective Value"] - \
                                                  gene_knockout_results["Knockout Objective Value"]
gene_knockout_results = gene_knockout_results.sort_values(by="Objective Value Change", ascending=False)

gene_knockout_results.to_csv("/mnt/NFS/fengch/TPM/drug/knock_Gal.csv", index=False)


