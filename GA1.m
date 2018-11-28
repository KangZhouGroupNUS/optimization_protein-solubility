%clear
rng default  % For reproducibility
fun = @fitness;
lb = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%options = optimoptions('ga','PlotFcn',@gaplotbestf,'PopulationSize',1000,'MaxGenerations',10000);
%rng default % For reproducibility
x = ga(fun,20,[],[],[],[],lb,[])