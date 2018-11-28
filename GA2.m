clear
rng default  % For reproducibility
%varnames = {'V1'; 'V2'; 'V3'; 'V4'; 'V5';'V6'; 'V7'; 'V8'; 'V9'; 'V10';'V11'; 'V12'; 'V13'; 'V14'; 'V15'; 'V16'; 'V17'; 'V18'; 'V19'; 'V20';'solubility'};
%Tbl = readtable('cleandata_1.csv','Filetype','text','ReadVariableNames',false);
%Tbl.Properties.VariableNames = varnames;
global CompactMdl;
CompactMdl=loadCompactModel('SVM');
%Mdl = fitrsvm(Tbl,'solubility','KernelFunction','gaussian','KernelScale','auto','Standardize',true)
fun = @fitness2;
%A = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1];
%b = [30;-30];
Aeq = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
beq = 30;
lb = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
ub = [30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30];
%opts = optimoptions('ga','PlotFcn',@gaplotbestf,'PopulationSize',1000,'MaxGenerations',10000);
%rng default % For reproducibility
%IntCon = 1;
%opts = optimoptions('ga','PlotFcn',@gaplotbestf,'PopulationSize',1000,'MaxGenerations',10000);
%x = ga(fun,20,A,b,[],[],lb,ub,[],IntCon,opts)
x = ga(fun,20,[],[],Aeq,beq,lb,ub)
output = fun(x)
sum = sum(x)


T = [];
rng default  % For reproducibility
global i;
global CompactMdl;
CompactMdl=loadCompactModel('SVM');
fun = @fitness2;
Aeq = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
beq = 30;
lb = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
ub = [30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30];
for i = 1 : 10
    x = ga(fun,20,[],[],Aeq,beq,lb,ub)
    output = fun(x)
    sum = sum(x)
    T = [T; i x output sum];
end
T