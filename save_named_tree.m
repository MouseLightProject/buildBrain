function save_named_tree(tree_mat_file_path, named_tree)
    save(tree_mat_file_path, '-mat', '-v7.3', 'named_tree') ;
end
