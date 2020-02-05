function result = is_component_a_binary_tree(A)
    % A should be the adjacency matrix for an undirected graph, and thus
    % symmetric.  All edges must also have the same weight.
    %
    % A must represent a *connected* graph, otherwise this test is not
    % reliable.
    
    if is_component_a_tree(A) ,
        degree = full(sum(A)) ;
        result = all(degree<=3) ;
    else
        result = false ;
    end
end
