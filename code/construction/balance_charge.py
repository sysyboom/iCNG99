from cobra.io import read_sbml_model
import pandas as pd

model = read_sbml_model('/mnt/NFS/fengch/new/models/merge40001_noc.xml')
charge_unbalanced_reactions_data = []
for reaction in model.reactions:
    reactant_total_charge = 0
    product_total_charge = 0
    for metabolite, stoichiometry in reaction.metabolites.items():
        if metabolite.charge is not None:
            if stoichiometry < 0:  
                reactant_total_charge += metabolite.charge * abs(stoichiometry)
            else:  
                product_total_charge += metabolite.charge * stoichiometry

    if reactant_total_charge != product_total_charge:
        charge_unbalanced_reactions_data.append({
            'Reaction ID': reaction.id,
            'Reactant Total Charge': reactant_total_charge,
            'Product Total Charge': product_total_charge
        })


charge_unbalanced_reactions_df = pd.DataFrame(charge_unbalanced_reactions_data)
with pd.ExcelWriter('/mnt/NFS/fengch/charge_merge40001_noc.xlsx') as writer:
    charge_unbalanced_reactions_df.to_excel(writer, sheet_name='Charge Unbalanced Reactions', index=False)


