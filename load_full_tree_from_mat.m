function [component_id, outtree] = load_full_tree_from_mat(full_tree_mat_file_path)
    load(full_tree_mat_file_path, '-mat', 'component_id', 'outtree') ;
    if ~exist('component_id', 'var') || ~exist('outtree', 'var') ,
        error('File %s is missing one or both of variables component_id and/or outtree', full_tree_mat_file_path) ;
    end
end
