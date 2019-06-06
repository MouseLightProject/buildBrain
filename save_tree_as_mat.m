function save_tree_as_mat(tree_mat_file_path, component_id, outtree)
    save(tree_mat_file_path, '-mat', '-v7.3', 'component_id', 'outtree') ;
end
