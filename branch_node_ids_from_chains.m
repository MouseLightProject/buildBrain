function result = branch_node_ids_from_chains(chains)
    chain_count = length(chains) ;
    heads_and_tails = zeros(chain_count, 2) ;
    for i = 1 : chain_count , 
        chain = chains{i} ;
        heads_and_tails(i,:) = [chain(1) chain(end)] ;
    end
    serial_heads_and_tails = heads_and_tails(:) ;
    nexus_node_ids = unique(serial_heads_and_tails) ;
    nexus_count = length(nexus_node_ids) ;
    degree_from_nexus_index = zeros(nexus_count, 1) ;
    for i = 1 : nexus_count ,
        node_id = nexus_node_ids(i) ;
        degree_from_nexus_index(i) = sum(serial_heads_and_tails==node_id) ;        
    end
    result = nexus_node_ids(degree_from_nexus_index==3) ;
end
