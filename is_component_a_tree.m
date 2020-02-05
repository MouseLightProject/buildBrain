function result = is_component_a_tree(A)
    % A should be the adjacency matrix for an undirected graph, and thus
    % symmetric.  All edges must also have the same weight.
    %
    % A must represent a *connected* graph, otherwise this test is not
    % reliable.
    
    node_count = size(A, 1) ;
    edge_count = full(sum(sum(A)))/2 ;
    result = (edge_count+1 == node_count) ;
end
