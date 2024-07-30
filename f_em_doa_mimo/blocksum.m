function L = blocksum(K,m,n)
% return sum of all block matrices of K with size m*n

[s1,s2] = size(K);
L = sum(reshape( sum( reshape(K,[],n,s2/n),3 ).', n, m, s1/m ),3).';

end