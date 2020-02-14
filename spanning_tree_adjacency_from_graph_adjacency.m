function [dA_tree, root] = spanning_tree_adjacency_from_graph_adjacency(A)
    % Turns a spanning tree of the undirected graph with adjacency matrix
    % A, as the adjacency matrix of the spanning tree.
    %
    % A should be the adjacency matrix for an undirected graph, and thus
    % symmetric.  All edges must also have unity weight.
    
    % Identify a leaf node to use as the root node, if possible.  Otherwise
    % use a min-degree node.
    degree = full(sum(A)) ;
    degree_modified = degree ;
    degree_modified(degree==0) = +inf ;
    [~, min_degree_node_index] = min(degree_modified) ;  % min degree is hopefully one, but sometimes not possible (e.g. if graph is just a loop)
    root = min_degree_node_index ;
    
    % Find the shortest path to the root from each node    
    node_count = size(A,1) ;
    [~,pred] = graphtraverse(A, root, 'Method', 'BFS', 'Directed', false) ;
    %pred = graphshortestpath_pred_only(A, root, 'directed', false, 'method', 'BFS') ;
    tree_ij_pairs = [((1:size(A,1))') (pred(:))] ;  % each row is a (child, parent)
    tree_ij_pairs(tree_ij_pairs(:,2)==0,:) = [];  % delete the row with the root as the child, and node 0 as the parent

    % Make the adjacency graph of the spanning tree
    dA_tree = sparse(tree_ij_pairs(:,1), tree_ij_pairs(:,2), 1, node_count, node_count) ; 
    %A_tree = max(dA_tree, dA_tree') ;
end
