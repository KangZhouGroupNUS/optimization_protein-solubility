%clear
nonlcon = @ellipsecons;
rng default  % For reproducibility

varnames = {'V1'; 'V2'; 'V3'; 'V4'; 'V5';'V6'; 'V7'; 'V8'; 'V9'; 'V10';'V11'; 'V12'; 'V13'; 'V14'; 'V15'; 'V16'; 'V17'; 'V18'; 'V19'; 'V20';'solubility'};
Tbl = readtable('cleandata_1.csv','Filetype','text','ReadVariableNames',false);
Tbl.Properties.VariableNames = varnames;
global CompactMdl;
CompactMdl=loadCompactModel('SVM');
%Mdl = fitrsvm(Tbl,'solubility','KernelFunction','gaussian','KernelScale','auto','Standardize',true)
fun = @svm_fitness;

Aeq = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
beq = 1;
lb = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
ub = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
options = optimoptions('ga','PlotFcn',@gaplotbestf,'PopulationSize',1000,'MaxGenerations',10000);
%rng default % For reproducibility
x = ga(fun,20,[],[],Aeq,beq,lb,ub,nonlcon,options)
output = fun(x)
[c, ceq]= ellipsecons(x)
sum = sum(x)
