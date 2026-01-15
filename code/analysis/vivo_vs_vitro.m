%% deal with model1(control group)
merge_after_vitro=importExcelModel('/mnt/NFS/fengch/new/models/paper/merge_after_YPD.xlsx')
model1=merge_after_vitro;
model1.rev=ones(numel(model1.rev),1);
model1.description = 'YPD';
model1=addUseVariablesDH(model1);
model1=addNetFluxVariablesNEW(model1);
% we want to add cellular growth which is at least 80 % of the maximum growth
sol=solveTFBAmodel(model1);
model1.var_lb(find(model1.f))=sol.val*0.8;
%% deal with model 2(experimental group)
merge_after_vivo=importExcelModel('/mnt/NFS/fengch/new/models/paper/merge_after_vivo.xlsx')
model2=merge_after_vivo;
model2.rev=ones(numel(model2.rev),1);
model2.description = 'vivo';
model2=addUseVariablesDH(model2);
model2=addNetFluxVariablesNEW(model2);
% we want to add cellular growth which is at least 80 % of the maximum growth%sol=solveTFBAmodel(model2);
sol=solveTFBAmodel(model2);
model2.var_lb(find(model2.f))=sol.val*0.8;

%% find up and downregulated reactions usig 5% up and down criteria

ratio_up=readmatrix('/mnt/NFS/fengch/TPM/vivo/Q3_up_value.xlsx');
ratio_down=readmatrix('/mnt/NFS/fengch/TPM/vivo/Q3_down_value.xlsx');
index_down=readmatrix('/mnt/NFS/fengch/TPM/vivo/Q3_down_index.xlsx');
index_up=readmatrix('/mnt/NFS/fengch/TPM/vivo/Q3_up_index.xlsx');
regRxnRatio=[ratio_up;ratio_down];
regIdx=[index_up;index_down];

%% create a Gex model which only intgerate relative expression
% add constraints for up and down regulated reactions
gex2=addRelExpCons2Models(model1,model2,model2.rxns(regIdx),regRxnRatio); % this can be also used addRelExpCons2Models
gex2.description='Regulated Model';
%sol=solveTFBAmodel_gurobi(gex2,false,'gurobi_direct','NA',false);
sol=solveTFBAmodel(gex2,false, 'cplex', false)
MCS=sol.val %% This is the maximum consistency score

%%  #### ALTERNATIVE ANALYSIS one should use findAltCombi function
% we shown in ecoli walk thorugh file
%time=500;
% findAltCombi(2,coM2,coM2.relExp.forB,path\_save,time);
path_save='/mnt/NFS/fengch/TPM/vivo/Q3.mat'
coM2=gex2;
coM2.var_lb(end)=MCS; % this will force maximum consistency can not be
%time=500;
% input argument : 1) numsol is number of alternatives
%                  2) model
%                  3) index for binary variables
%                  4) path for saving the resultes
%                  5) time for solver
%findAltCombi(1000,coM2,coM2.relExp.forB,path\_save);
%     f=zeros(numel(model.f),1);
%     f(index)=1;
%     model.f=f;
index=[coM2.objIndex1B]
%findAltCombi_Gurobi(100,coM2,index,path\_save,43200)
findAltCombi(100,coM2,index,path_save)
