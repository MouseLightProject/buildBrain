function result = branch_node_ids_from_edges(edges)
    serial_edges = edges(:) ;
    node_ids = unique(serial_edges) ;
    node_count = length(node_ids) ;
    degree_from_node_index = zeros(node_count, 1) ;
    for i = 1 : node_count ,
        node_id = node_ids(i) ;
        degree_from_node_index(i) = sum(serial_edges==node_id) ;        
    end
    result = node_ids(degree_from_node_index==3) ;
end
