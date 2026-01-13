import cobra
import pytest


# 定义 model fixture
@pytest.fixture
def model():
    # 加载模型
    model = cobra.io.read_sbml_model('/mnt/NFS/fengch/new/models/merge_after.xml')
    return model


# 测试敲除基因是否会影响目标函数
def test_gene_knockout(model):
    gene_id_to_knockout = 'CNAG_06508'

    # 检查基因是否存在
    try:
        gene = model.genes.get_by_id(gene_id_to_knockout)
    except KeyError:
        print(f"Gene {gene_id_to_knockout} not found in the model.")
        return False

    # 敲除基因并优化模型
    model.genes.get_by_id(gene_id_to_knockout).knock_out()
    solution = model.optimize()

    # 输出目标函数流值
    print(f"Objective value after knocking out {gene_id_to_knockout}: {solution.objective_value}")

    # 判断是否有流
    if solution.objective_value > 1e-6:
        print(f"After knocking out {gene_id_to_knockout}, the network still has flux through the objective function.")
        return True
    else:
        print(f"After knocking out {gene_id_to_knockout}, there is no flux through the objective function.")
        return False
