function result = asssign_random_colors(named_trees)
    tree_count = length(named_trees) ;
    color_map = distinct_hues_simple() ;    
    rep_count = ceil(tree_count/size(color_map,1)) ;
    color_from_tree_index = repmat(color_map, [rep_count 1]) ;
    result = named_trees ;
    for tree_index = 1 : tree_count ,
        result(tree_index).color = color_from_tree_index(tree_index,:) ;
    end
end
