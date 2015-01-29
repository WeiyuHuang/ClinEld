function [ ec ] = eigen_centrality( A )
%compute closeness centrality given an adjacency matrix A

[ec, ~] = eigs(A, 1);

if sum(ec) < 0
    ec = -ec;
end

ec = ec';

end

