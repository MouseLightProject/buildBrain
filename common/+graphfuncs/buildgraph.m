function tree_ij_pairs = buildgraph(A, root)
    % Turns a graph with an indicated root node into a rooted tree, rooted
    % at the given root.
    % In the output, the path from a node to the root will be the shortest
    % such path in the orginal graph.
    % 
    % Asub should be the adjacency matrix for an undirected graph, and thus
    % symmetric.
    %
    % The output is a (node_count-1) x 2 matrix.
    % Each row lists a (child, parent) pair.
    % The root node never appears in the 1st column, so the output has
    % node_count-1 rows.
    
%     if nargin<2 ,
%         source = find(sum(Asub)==1, 1) ;
%     end
    [~,~,pred] = graphshortestpath(A, root, 'directed', false) ;
    tree_ij_pairs = [((1:size(A,1))') (pred(:))] ;  % each row is a (child, parent)
    tree_ij_pairs(tree_ij_pairs(:,2)==0,:) = [];  % delete the row with the root as the child, and node 0 as the parent
end
