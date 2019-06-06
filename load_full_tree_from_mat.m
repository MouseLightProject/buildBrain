function [component_id, outtree] = load_full_tree_from_mat(full_tree_mat_file_path)
    load(full_tree_mat_file_path, '-mat', 'component_id', 'outtree') ;
end
