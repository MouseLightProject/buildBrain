function is_pool_tree_near_query_tree_from_pool_tree_id = ...
        find_trees_near_query_tree(query_tree_as_named_tree, ...
                                   pool_kd_tree, ...
                                   pool_node_count, ...
                                   pool_tree_id_from_pool_node_id, ...
                                   pool_tree_count, ...
                                   blast_radius_in_um)

    query_tree_node_count = size(query_tree_as_named_tree.xyz, 1) ;
                               
    xyz_from_query_tree_node_id = query_tree_as_named_tree.xyz ;
    
    nearby_pool_node_ids_from_query_tree_node_id = rangesearch(pool_kd_tree, xyz_from_query_tree_node_id, blast_radius_in_um) ;
    
    is_pool_tree_node_near_query_tree_from_pool_node_id = false(pool_node_count, 1) ;
    for query_tree_node_id = 1 : query_tree_node_count ,
        nearby_pool_node_ids = nearby_pool_node_ids_from_query_tree_node_id{query_tree_node_id} ;
        is_pool_tree_node_near_query_tree_from_pool_node_id(nearby_pool_node_ids) = true ;
    end
    pool_tree_ids_near_query_tree = pool_tree_id_from_pool_node_id(is_pool_tree_node_near_query_tree_from_pool_node_id) ;  % may contain repeats
    
    is_pool_tree_near_query_tree_from_pool_tree_id = false(pool_tree_count, 1) ;
    is_pool_tree_near_query_tree_from_pool_tree_id(pool_tree_ids_near_query_tree) = true ;
end
