close all
clear

sample_date = '2019-05-27' ;
input_swcs_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/shared_tracing/Finished_Neurons/%s', sample_date) ;
input_sample_folder_path = sprintf('/nrs/mouselight/SAMPLES/%s', sample_date) ;
input_trees_folder_path = ...
    sprintf('/nrs/mouselight/cluster/classifierOutputs/%s/build-brain-output/full-as-mats', sample_date) ;
output_swcs_folder_path = ...
    sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/build-brain-output/full-trees-near-interesting-somata', sample_date) ;

% Get the metadata for the rendered sample
[shape_xyz, origin, spacing, jaws_origin] = load_sample_shape_origin_and_spacing(input_sample_folder_path);

% Load the fragments from disk                                                
tic_id = tic() ;
full_trees_as_named_trees_in_erhan_coords = load_forest_from_folder_of_mats(input_trees_folder_path) ;
elapsed_time = toc(tic_id) ;
fprintf('Elapsed time to load in all trees was %g seconds\n', elapsed_time) ;

% Extract all the fragment nodes, along with a look-up table to map them back to
% the fragment they came from
full_tree_count = length(full_trees_as_named_trees_in_erhan_coords) ;
[xyz_in_erhan_coords_from_overall_node_id, tree_id_from_overall_node_id] = ...
    point_cloud_from_forest_as_named_trees(full_trees_as_named_trees_in_erhan_coords) ;
tree_node_count = length(tree_id_from_overall_node_id) ;

% Check that the fragment nodes are on the grid as we expect
%     all_fragment_ijks_unrounded = (all_fragment_xyzs_in_erhan_coords - origin) ./ spacing + 0.5 ;  % these should be one-based ijks
%     all_fragment_ijks = round(all_fragment_ijks_unrounded) ;  % these should be one-based ijks
%     fragment_offsets = all_fragment_ijks_unrounded - all_fragment_ijks ;
%     assert( all( all( abs(fragment_offsets) < spacing/1024 ) ) ) ;  % the offsets should just be floating-point error                              
fraction_of_uncorrected_fragment_nodes_at_voxel_centers = ...
    fraction_of_xyzs_at_voxel_centers_using_erhan_conventions(xyz_in_erhan_coords_from_overall_node_id, spacing, origin)
assert( fraction_of_uncorrected_fragment_nodes_at_voxel_centers == 1) ;

% Convert fragment node xyz coords to voxel indices
unrounded_ijk1_from_overall_node_id = (xyz_in_erhan_coords_from_overall_node_id - origin) ./ spacing + 0.5 ;  % these should be one-based ijks
ijk1_from_overall_node_id = round(unrounded_ijk1_from_overall_node_id) ;  % these should be one-based ijks

% % Check the alignment of the fragments to the fluorescence signal
% desired_sample_count = 100 ;
% [f,a] = figure_to_check_alignment_to_stack(sample_folder_path, ijk1_from_fragment_node_id, desired_sample_count) ;

% Convert one-based voxel indices to voxel center coords in JaWS coord system
% all_uncorrected_fragment_xyzs = jaws_origin + spacing .* (all_uncorrected_fragment_ijks-1) ;  % um, n x 3
% all_fragment_xyzs = jaws_origin + spacing .* (all_fragment_ijks-1) ;  % um, n x 3
xyz_from_overall_node_id = jaws_origin + spacing .* (ijk1_from_overall_node_id-1) ;  % um, n x 3
    % ijk1_from_fragment_node_id are 1-based, we want to map them to voxel centers
    % in JaWS coordinates

% Load in the "neurons" specified by Jayaram
somata_count = 5 ;
somata_as_named_trees = empty_named_tree_struct(somata_count,1) ;
somata_as_named_trees(1).name = 'G-091-soma' ;
somata_as_named_trees(1).xyz = [73402.0, 17144.8, 37251.4] ;
somata_as_named_trees(2).name = 'G-082-soma' ;
somata_as_named_trees(2).xyz = [74814.0, 15426.1, 38458.4] ;
somata_as_named_trees(3).name = 'G-074-soma' ;
somata_as_named_trees(3).xyz = [73482.8, 16474.8, 35215.7] ;
somata_as_named_trees(4).name = 'G-031-soma' ;
somata_as_named_trees(4).xyz = [78269.0, 16192.7, 34378.1] ;
somata_as_named_trees(5).name = 'G-008-soma' ;
somata_as_named_trees(5).xyz = [74373.7, 15286.3, 37962.8] ;
for somata_index = 1 : somata_count ,
    somata_as_named_trees(somata_index).color = [1 1 1] ;
    somata_as_named_trees(somata_index).r = 1 ;
    somata_as_named_trees(somata_index).tag_code = 0 ;
    somata_as_named_trees(somata_index).parent = -1 ;
end

forest_kd_tree  = KDTreeSearcher(xyz_from_overall_node_id) ;        
blast_radius_in_um = 20 ;  % um
is_tree_near_somata_from_tree_id = ...
    find_trees_near_query_trees(...
        somata_as_named_trees, ...
        forest_kd_tree, ...
        tree_node_count, ...
        tree_id_from_overall_node_id, ...
        full_tree_count, ...
        blast_radius_in_um) ;

% Filter the fragments, convert to jaws coords
nearby_full_trees_in_erhan_coords = full_trees_as_named_trees_in_erhan_coords(is_tree_near_somata_from_tree_id) ;
nearby_full_trees = convert_named_trees_to_jaws_coords(nearby_full_trees_in_erhan_coords, origin, spacing, jaws_origin) ;
  
%%
% Write the swc's that don't exist already                          
save_swcs_from_named_trees(output_swcs_folder_path, ...                         
                           nearby_full_trees) ;
