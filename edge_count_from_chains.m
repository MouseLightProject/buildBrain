function result = edge_count_from_chains(chains)
    %chain_count = length(chains) ;
    node_count_from_chain_index = cellfun(@length, chains) ;
    edge_count_from_chain_index = node_count_from_chain_index-1 ;
    result = sum(edge_count_from_chain_index) ;
end
