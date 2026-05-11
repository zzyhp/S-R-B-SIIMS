function [idx, C] = kmeans1(T, K)
T = double(T);

X = T(:);     

[idx, C] = kmeans(X, K);

end