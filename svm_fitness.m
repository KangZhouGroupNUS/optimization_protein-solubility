function y = svm_fitness(x)
global CompactMdl;
y = -predict(CompactMdl,x);
end

