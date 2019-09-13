function y = fitness5(x)
global CompactMdl;
global i;
Tb2 = readtable('matlab_5.csv','Filetype','text','ReadVariableNames',false);
initial_sample= Tb2(i,1:20);
initial=table2array(initial_sample);
%length_table=[266 271 289 275 133 303 155 270 128 141];
length_table=[509 543 388 459 620 589];
length = length_table(i);
x1=(initial*length+x)/(20+length);
y = -predict(CompactMdl,x1);
end

