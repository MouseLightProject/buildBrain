close all
clear

input_folder_path = '/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2019-05-27/build-brain-output/frags-as-mats' ;
output_file_path = '/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2019-05-27/build-brain-output/frags.mat' ;
maximum_core_count_desired = inf ;

% Get the pool ready
poolobj = gcp('nocreate');  % If no pool, do not create new one.
if isempty(poolobj) ,
    parpool([1 maximum_core_count_desired]) ;
end
poolobj = gcp('nocreate');  % If no pool, do not create new one.
core_count = poolobj.NumWorkers ;
fprintf('Using %d cores.\n', core_count) ;

frag_mat_file_names = simple_dir(fullfile(input_folder_path, 'auto-cc-*-fragments.mat')) ;
tree_count = length(frag_mat_file_names) ;
fragments_from_tree_id = cell(1, tree_count) ;
tic_id = tic() ;
parfor_progress(tree_count) ;
parfor tree_index = 1 : tree_count ,
    frag_mat_file_name = frag_mat_file_names{tree_index} ;
    frag_mat_file_path = fullfile(input_folder_path, frag_mat_file_name) ;
    fragments_from_tree_id{tree_index} = load_anonymous(frag_mat_file_path) ;
    parfor_progress() ;
end
parfor_progress(0) ;
elapsed_time = toc(tic_id) ;
fprintf('Elapsed time to read in all fragment .mats was %g seconds\n', elapsed_time) ;

fragments_as_swc_arrays = horzcat(fragments_from_tree_id{:}) ;  

tic_id = tic() ;
save(output_file_path, 'fragments_as_swc_arrays', '-v7.3', '-nocompression') ;
elapsed_time = toc(tic_id) ;
fprintf('Elapsed time to save in all fragments in single .mat was %g seconds\n', elapsed_time) ;
