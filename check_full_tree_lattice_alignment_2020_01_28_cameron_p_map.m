sample_date = '2020-01-28-cameron-p-map' ;
input_mats_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/build-brain-output/full-as-named-tree-mats', sample_date) ;
input_sample_folder_path = '/nrs/mouselight-v/lightsheet/2020-01-28' ;
input_fragments_file_path = ...
    sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/build-brain-output/frags.mat', sample_date) ;

% Get the metadata for the rendered sample
[shape_xyz, origin, spacing, jaws_origin] = load_sample_shape_origin_and_spacing(input_sample_folder_path);

% Load the fragments from disk                                                
tic_id = tic() ;
forest = load_forest_from_folder_of_named_tree_mats(input_mats_folder_path) ;
elapsed_time = toc(tic_id) ;
fprintf('Elapsed time to load in all fragments from a single .mat was %g seconds\n', elapsed_time) ;

% Extract all the fragment nodes, along with a look-up table to map them back to
% the fragment they came from
tree_count = length(forest) 
[xyz_in_erhan_coords_from_tree_node_id, fragment_id_from_tree_node_id] = point_cloud_from_forest_as_named_trees(forest) ;
tree_node_count = length(fragment_id_from_tree_node_id) 

% Check that the fragment nodes are on the grid as we expect
%     all_fragment_ijks_unrounded = (all_fragment_xyzs_in_erhan_coords - origin) ./ spacing + 0.5 ;  % these should be one-based ijks
%     all_fragment_ijks = round(all_fragment_ijks_unrounded) ;  % these should be one-based ijks
%     fragment_offsets = all_fragment_ijks_unrounded - all_fragment_ijks ;
%     assert( all( all( abs(fragment_offsets) < spacing/1024 ) ) ) ;  % the offsets should just be floating-point error                              
fraction_of_uncorrected_fragment_nodes_at_voxel_centers = ...
    fraction_of_xyzs_at_voxel_centers_using_erhan_conventions(xyz_in_erhan_coords_from_tree_node_id, spacing, origin)
assert( fraction_of_uncorrected_fragment_nodes_at_voxel_centers == 1) ;

% Convert fragment node xyz coords to voxel indices
unrounded_ijk1_from_tree_node_id = (xyz_in_erhan_coords_from_tree_node_id - origin) ./ spacing + 0.5 ;  % these should be one-based ijks
ijk1_from_tree_node_id = round(unrounded_ijk1_from_tree_node_id) ;  % these should be one-based ijks

% Check the alignment of the fragments to the fluorescence signal
desired_sample_count = 100 ;
[f,a] = figure_to_check_alignment_to_stack(input_sample_folder_path, ijk1_from_tree_node_id, desired_sample_count) ;
