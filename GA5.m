T = [];
rng default  % For reproducibility
global i;
global CompactMdl;
CompactMdl=loadCompactModel('SVM6');
fun = @fitness5;
Aeq = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
beq = 20;
% amino acid sequence A R N D C E Q G H I L K M F P S T W Y V
lb = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
ub = [20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20];
for i = 1 : 6
    x = ga(fun,20,[],[],Aeq,beq,lb,ub)
    output = fun(x)
    sum = sum(x)
    T = [T; i x output sum];
end
T