function L = convsum(K,m,n)
% K is devided into m*n block matrices, then sum over elements in each
% block

[s1,s2] = size(K);
L = reshape( permute( reshape( K, [], n, s2/n ), [2,1,3] ), [n,m,s1/m,s2/n]);
L = sum(sum(L,2),1);
L = reshape(L,s1/m,s2/n);



end