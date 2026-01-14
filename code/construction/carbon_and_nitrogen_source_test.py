import cobra
import pandas as pd

# Alternative nitrogen sources
exchange_reactions_group1 = []
exchange_reactions_group2 = ['EX_040', 'EX_111', 'EX_112', 'EX_113', 'EX_114','EX_115','EX_116','EX_117','EX_118','EX_119','EX_120','EX_121','EX_122','EX_123','EX_124','EX_125','EX_126','EX_127','EX_128','EX_129','EX_130']
fixed_reaction_id = 'EX_015' 

biomass1_reaction_id = 'biomass2_c' # with capsule
biomass2_reaction_id = 'biomass1_c' # no capsule

results = []

for i in range(len(exchange_reactions_group1)):
    print(f"Starting iteration {i + 1} for Group 1 (biomass1)...")

    model = cobra.io.read_sbml_model('/mnt/NFS/fengch/new/models/mergemodeltest40001.xml')
    model.solver = 'cplex'

    current_reaction_id = exchange_reactions_group1[i]
    
    reaction = model.reactions.get_by_id(current_reaction_id)
    reaction.lower_bound = -1000
    reaction.upper_bound = 1000

    fixed_reaction = model.reactions.get_by_id(fixed_reaction_id)
    fixed_reaction.lower_bound = 0
    fixed_reaction.upper_bound = 0

    target_reaction = model.reactions.get_by_id(biomass1_reaction_id)
    model.objective = target_reaction
    model.objective_direction = 'max'

    solution = model.optimize()

    results.append({
        'Iteration': i + 1,
        'Group': 'Group 1',
        'Modified_Reaction': current_reaction_id,
        'Fixed_Reaction': fixed_reaction_id,
        'Objective_Function': biomass1_reaction_id,
        'Optimized_Value': solution.objective_value
    })

    print(f'Iteration {i + 1} - {current_reaction_id} and {fixed_reaction_id} modified:')
    print('Optimized value:', solution.objective_value)

for i in range(len(exchange_reactions_group2)):
    print(f"Starting iteration {i + 1} for Group 2 (biomass2)...")

    model = cobra.io.read_sbml_model('/mnt/NFS/fengch/new/models/mergemodeltest40001.xml')
    model.solver = 'cplex'

    current_reaction_id = exchange_reactions_group2[i]

    reaction = model.reactions.get_by_id(current_reaction_id)
    reaction.lower_bound = -1000
    reaction.upper_bound = 1000

    fixed_reaction = model.reactions.get_by_id(fixed_reaction_id)
    fixed_reaction.lower_bound = 0
    fixed_reaction.upper_bound = 0

    target_reaction = model.reactions.get_by_id(biomass2_reaction_id)
    model.objective = target_reaction
    model.objective_direction = 'max'

    solution = model.optimize()

    results.append({
        'Iteration': i + 1,
        'Group': 'Group 2',
        'Modified_Reaction': current_reaction_id,
        'Fixed_Reaction': fixed_reaction_id,
        'Objective_Function': biomass2_reaction_id,
        'Optimized_Value': solution.objective_value
    })

    print(f'Iteration {i + 1} - {current_reaction_id} and {fixed_reaction_id} modified:')
    print('Optimized value:', solution.objective_value)

df_results = pd.DataFrame(results)

df_results.to_excel('/mnt/NFS/fengch/new/transport/nitrogen_biomass_results.xlsx', index=False)

print("Results have been saved to 'results.xlsx'.")
