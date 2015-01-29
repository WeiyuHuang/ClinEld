function [ N ] = normalized_adjacencyMatrix( A )
%Given A as an adjacency matrix, compute N as the normalized adjacency
%matrix where entry (i,j) of N is defined as A(i,j) / sqrt(d(i) * d(j)) and
%d(i) is the degree of A, i.e. d(i) = \sum_j A(i,j)

D = sum(A);
N = zeros(size(A));
for i = 1:size(A, 1);
    for j = i:size(A, 1);
        N(i,j) = A(i,j) ./ sqrt(D(i) * D(j));
        N(j,i) = N(i,j);
    end
end


end

