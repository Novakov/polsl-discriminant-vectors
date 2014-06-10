function [d,gamma] = runFile(fileName, K)

dataset = csvread(fileName);
x = {};

x{1} = dataset(dataset(:,end)==0, 1:end-1);
x{2} = dataset(dataset(:,end)==1, 1:end-1);

[d, gamma] = discrimVectors(x, K, 0.5);

bar(abs(transpose(d)))