function y = fitness2(x)
global CompactMdl;
global i;
Tb2 = readtable('sample_10.csv','Filetype','text','ReadVariableNames',false);
initial_sample= Tb2(i,1:20);
initial=table2array(initial_sample);
%initial = [0.08285714 0.07142857 0.02285714 0.03714286 0.02000000 0.06571429 0.05714286 0.11142857 0.01714286 0.04571429 0.12000000 0.02285714 0.02857143 0.03142857 0.06285714 0.06285714 0.04000000 0.02285714 0.02285714 0.05428571]; 
length_table=[266 271 289 275 133 280 155 270 296 160];
length = length_table(i);
x1=(initial*length+x)/(20+length);
y = -predict(CompactMdl,x1);
end

