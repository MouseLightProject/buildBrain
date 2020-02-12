function node_ids_from_chain_id = chains_from_tree(A)
    % Given symmetric adjacency matrix A representing an undirected tree,
    % returns a set of chains.  The first and last node of each chain is a
    % node_id into A, and each is a nexus node.  All nodes in A are
    % included in at least one chain.

%     % For debugging
%     G = graph(dA) ;
%     figure() ;
%     plot(G) ;

    if isempty(A) ,
        node_ids_from_chain_id = cell(1,0) ;
        return
    end

    % Get a rooted directed tree from the undirected tree
    % The 'root' node will basically be a random leaf, which is what we
    % want.
    dA = spanning_tree_adjacency_from_graph_adjacency(A) ;
    
    % Find the non-root nexus nodes---each of these will be the start of a
    % chain
    in_degree_from_node_id = sum(dA,1)' ;
    out_degree_from_node_id = sum(dA,2) ;
    is_root_from_node_id = (out_degree_from_node_id==0) ;
    is_chain_node_from_node_id = (in_degree_from_node_id==1) & (out_degree_from_node_id==1) ;
    is_nexus_from_node_id = ~is_chain_node_from_node_id ;
    is_nonroot_nexus_from_node_id = is_nexus_from_node_id & ~is_root_from_node_id ;
    distal_nexus_node_id_from_chain_id = find(is_nonroot_nexus_from_node_id) ;
    chain_count = length(distal_nexus_node_id_from_chain_id) ;
    
    % Trace each chain from its distal nexus node to its proximal nexus
    % node
    node_ids_from_chain_id = cell(chain_count, 1) ;
    for chain_id = 1 : chain_count ,
        distal_node_id = distal_nexus_node_id_from_chain_id(chain_id) ;
        [node_ids_in_this_chain, proximal_nexus_node_id] = ...
            trace_chain_in_directed(dA, distal_node_id) ;
        node_ids_from_chain_id{chain_id} = [node_ids_in_this_chain proximal_nexus_node_id] ;
    end
end
