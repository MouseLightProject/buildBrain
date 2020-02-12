function edges = edges_from_chains(node_ids_from_chain_id)
    chain_count = length(node_ids_from_chain_id) ;
    edges_from_chain_id = cell(chain_count, 1) ;
    for chain_id = 1 : chain_count ,
        node_ids = node_ids_from_chain_id{chain_id} ;
        to = node_ids(1:end-1) ;
        from = node_ids(2:end) ;
        edges_from_chain_id{chain_id} = [from(:) to(:)]' ;  % 2 x edge_count_for_this_chain
    end
    edges = [edges_from_chain_id{:}]' ;  % edge_count x 2, all the edges, across all chains
end
