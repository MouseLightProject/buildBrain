function [A_decimated, xyz_decimated, original_node_id_from_decimated_node_id] = decimate_tree(A_tree, xyz, desired_sampling_interval)
    % Decimates a tree represented as an undirected adjacency matrix
    % A_tree, with nodes at coords given by xyz, in um.  desired_sampling_interval
    % should be in um.
    %
    % Result is returned as an undirected tree graph, A_decimated, at roughly the desired
    % sampling interval.  The xyz coords of the nodes of A_decimated are in
    % the rows of xyz_decimated.  Each node of A_decimated is a node of
    % A_tree, and original_node_id_from_decimated_node_id allows the called
    % to map between them.
    
    node_count = size(A_tree,1) ;
    if node_count==0 ,
        % Special case for empty tree
        A_decimated = A_tree ;
        xyz_decimated = xyz ;
        original_node_id_from_decimated_node_id = [] ;
    else
        % The typical, non-annoying case
        
        % Decompose tree into chains
        node_ids_from_chain_id = chains_from_tree(A_tree) ;

        % Decimate chains
        decimated_node_ids_from_chain_id = decimate_chains(node_ids_from_chain_id, desired_sampling_interval) ;

        % Convert chains to edges
        edges_using_decimated_node_ids = edges_from_chains(decimated_node_ids_from_chain_id) ;

        % Defragment the node ids
        [edges_using_new_node_ids, xyz_decimated, original_node_id_from_decimated_node_id] = ...
            defragment_node_ids_in_edges(edges_using_decimated_node_ids, xyz) ;

        % Convert edges to (sparse) adjacency
        decimated_node_count = length(original_node_id_from_decimated_node_id) ;
        A_decimated = undirected_adjacency_from_edges(edges_using_new_node_ids, decimated_node_count) ;    
    end
end
