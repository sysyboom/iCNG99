import cobra
import pytest

def model():
    model = cobra.io.read_sbml_model('/mnt/NFS/fengch/new/models/merge_after.xml')
    return model

def test_gene_knockout(model):
    gene_id_to_knockout = 'CNAG_06508'

    try:
        gene = model.genes.get_by_id(gene_id_to_knockout)
    except KeyError:
        print(f"Gene {gene_id_to_knockout} not found in the model.")
        return False

    model.genes.get_by_id(gene_id_to_knockout).knock_out()
    solution = model.optimize()

    print(f"Objective value after knocking out {gene_id_to_knockout}: {solution.objective_value}")

    if solution.objective_value > 1e-6:
        print(f"After knocking out {gene_id_to_knockout}, the network still has flux through the objective function.")
        return True
    else:
        print(f"After knocking out {gene_id_to_knockout}, there is no flux through the objective function.")
        return False
