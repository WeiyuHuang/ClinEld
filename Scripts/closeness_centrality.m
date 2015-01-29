function [ cc ] = closeness_centrality( A )
%compute closeness centrality given an adjacency matrix A

D = ones(size(A)) ./ A;
graph_shortestPath = graphallshortestpaths(sparse(D));
graph_shortestPath = graph_shortestPath - diag(diag(graph_shortestPath));

cc_reci = sum(graph_shortestPath);
cc = ones(size(cc_reci))./cc_reci;

end

