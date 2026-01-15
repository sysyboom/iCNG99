%%% heat data 1
%% deal with model1(control group)
merge_after_YPD_heat=importExcelModel('/mnt/NFS/fengch/new/models/paper/merge_after_YPD_heat.xlsx')
model1=merge_after_YPD_heat; % gal
model1.rev=ones(numel(model1.rev),1);
model1.description = 'gal';
model1=addUseVariablesDH(model1);
model1=addNetFluxVariablesNEW(model1);
% we want to add cellular growth which is at least 80 % of the maximum growth
sol=solveTFBAmodel(model1);
model1.var_lb(find(model1.f))=sol.val*0.8;
%% deal with model 2(experimental group)
model2=merge_after_YPD_heat; % glu
model2.rev=ones(numel(model2.rev),1);
model2.description = 'gal';
model2=addUseVariablesDH(model2);
model2=addNetFluxVariablesNEW(model2);
% we want to add cellular growth which is at least 80 % of the maximum growth%sol=solveTFBAmodel(model2);
sol=solveTFBAmodel(model2);
model2.var_lb(find(model2.f))=sol.val*0.8;

%% find up and downregulated reactions usig 5% up and down criteria

ratio_up=readmatrix('/mnt/NFS/fengch/TPM/heat/toler2_up_value.xlsx');
ratio_down=readmatrix('/mnt/NFS/fengch/TPM/heat/toler2_down_value.xlsx');
index_down=readmatrix('/mnt/NFS/fengch/TPM/heat/toler2_down_index.xlsx');
index_up=readmatrix('/mnt/NFS/fengch/TPM/heat/toler2_up_index.xlsx');
regRxnRatio=[ratio_up;ratio_down];
regIdx=[index_up;index_down];

%% create a Gex model which only intgerate relative expression
% add constraints for up and down regulated reactions
gex2=addRelExpCons2Models(model1,model2,model2.rxns(regIdx),regRxnRatio); % this can be also used addRelExpCons2Models
gex2.description='Regulated Model';
%%
regIdx_met=readmatrix('/mnt/NFS/fengch/new_data/new_results_YPD/index2.xlsx')
regMetRatio=readmatrix('/mnt/NFS/fengch/new_data/new_results_YPD/meta2.xlsx')
gex2.S = model2.S;
gex2.mets = model2.mets;
gex2.rxns = model2.rxns;
tmpModel=addRelMetabolite_coef(gex2,gex2,model2.mets(regIdx_met),regMetRatio,true)
sol=solveTFBAmodel(tmpModel);

MCS=sol.val %% This is the maximum consistency score


%%  #### ALTERNATIVE ANALYSIS one should use findAltCombi function
% we shown in ecoli walk thorugh file
%time=500;
% findAltCombi(2,coM2,coM2.relExp.forB,path_save,time);
path_save='/mnt/NFS/fengch/TPM/heat/heat_YPD_2.mat'
coM2 = tmpModel;
coM2.var_lb(end)=MCS; % this will force maximum consistency can not be 
%time=500;
% input argument : 1) numsol is number of alternatives
%                  2) model
%                  3) index for binary variables 
%                  4) path for saving the resultes
%                  5) time for solver
%findAltCombi(1000,coM2,coM2.relExp.forB,path_save);
%     f=zeros(numel(model.f),1);
%     f(index)=1;
%     model.f=f;
index=[coM2.objIndex1B;coM2.metB];
findAltCombi_Gurobi(100,coM2,index,path_save,43200)

% save('iShiiModelThermoRelative.mat','storeModels')

