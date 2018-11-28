varnames = {'V1'; 'V2'; 'V3'; 'V4'; 'V5';'V6'; 'V7'; 'V8'; 'V9'; 'V10';'V11'; 'V12'; 'V13'; 'V14'; 'V15'; 'V16'; 'V17'; 'V18'; 'V19'; 'V20';'solubility'};
Tbl = readtable('cleandata_1.csv','Filetype','text','ReadVariableNames',false);
Tbl.Properties.VariableNames = varnames;
SVM = fitrsvm(Tbl,'solubility','OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus'))
%Cost (in R) refers to Box Constraint (MATLAB) 1.26
%1/sqrt(Gamma) refers to Kernel Scale (MATLAB) 5.8722
%Epsilon 0.057
%For the optimization and plots, the objective function is log(1 + cross-validation loss) for regression, and the misclassification rate for classification.
saveCompactModel(SVM,'SVM');
%CompactMdl = loadCompactModel('SVM');