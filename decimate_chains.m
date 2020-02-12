function decimated_node_ids_from_chain_id = decimate_chains(node_ids_from_chain_id, xyz, desired_sampling_interval)
    % Decimate chains
    chain_count = length(node_ids_from_chain_id) ;
    decimated_node_ids_from_chain_id = cell(chain_count, 1) ;
    for chain_id = 1 : chain_count ,
        % Get the node ids for this chain
        node_ids_in_chain = node_ids_from_chain_id{chain_id} ;               
        decimated_node_ids_in_chain = decimate_chain(node_ids_in_chain, xyz, desired_sampling_interval) ;
        decimated_node_ids_from_chain_id{chain_id} = decimated_node_ids_in_chain ;
    end    
end
