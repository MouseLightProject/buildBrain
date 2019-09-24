function convert_full_trees_to_swc(trees_as_swc_folder_path, trees_as_mat_folder_path)
    tic_id = tic() ;
    if ~exist(trees_as_swc_folder_path, 'file') ,
        mkdir(trees_as_swc_folder_path) ;
    end
    mat_file_names = simple_dir(fullfile(trees_as_mat_folder_path, 'auto*.mat')) ;
    file_count = length(mat_file_names) ;
    
    %raw_color_from_file_index = jet(file_count);
    %color_from_file_index = raw_color_from_file_index(randperm(file_count),:);
    
    color_map = distinct_hues_simple() ;    
    rep_count = ceil(file_count/size(color_map,1)) ;
    color_from_file_index = repmat(color_map, [rep_count 1]) ;
    
    fprintf('Converting %d trees as .mat to .swc...\n', file_count) ;
    parfor_progress(file_count) ;
    for i = 1 : file_count ,
        mat_file_name = mat_file_names{i} ;
        swc_base_name = base_name_from_file_name(mat_file_name) ;
        swc_file_name = horzcat(swc_base_name, '.swc') ;
        swc_file_path = fullfile(trees_as_swc_folder_path, swc_file_name) ;
        if ~exist(swc_file_path, 'file') ,
            mat_file_path = fullfile(trees_as_mat_folder_path, mat_file_name) ;
            mat_contents = load('-mat', mat_file_path) ;
            tree_as_struct = mat_contents.outtree ;
            %component_id = mat_contents.component_id ;
            tree_as_swc_array = tree_as_swc_array_from_tree_as_struct(tree_as_struct) ;
            color = color_from_file_index(i, :) ;
            neuron_name = swc_base_name ;
            save_swc(swc_file_path, tree_as_swc_array, neuron_name, color) ;
        end
        parfor_progress() ;
    end
    parfor_progress(0) ;
    elapsed_time = toc(tic_id) ;
    fprintf('Elapsed time for convert_full_trees_to_swc() was %g seconds\n', elapsed_time) ;
end
