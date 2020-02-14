function tree_ij_pairs = rooted_spanning_tree_edges_from_graph(A, root)
    % Turns a graph with an indicated root node into a rooted tree, rooted
    % at the given root.
    % In the output, the path from a node to the root will be the shortest
    % such path in the orginal graph.
    % 
    % A should be the adjacency matrix for an undirected graph, and thus
    % symmetric.  All edges must also have the same weight.
    %
    % The output is a (node_count-1) x 2 matrix.
    % Each row lists a (child, parent) pair.
    % The root node never appears in the 1st column, so the output has
    % node_count-1 rows.
    
    if nargin<2 ,
        degree = full(sum(A)) ;
        degree_modified = degree ;
        degree_modified(degree==0) = +inf ;
        [~, min_degree_index] = min(degree_modified) ;  % min degree is hopefully one, but sometimes not possible (e.g. if graph is just a loop)
        root = min_degree_index ;
    end
    %pred = graphshortestpath_pred_only(A, root, 'directed', false, 'method', 'BFS') ;
    [~,pred] = graphtraverse(A, root, 'Method', 'BFS', 'Directed', false) ;
    tree_ij_pairs = [((1:size(A,1))') (pred(:))] ;  % each row is a (child, parent)
    tree_ij_pairs(tree_ij_pairs(:,2)==0,:) = [];  % delete the row with the root as the child, and node 0 as the parent
end
