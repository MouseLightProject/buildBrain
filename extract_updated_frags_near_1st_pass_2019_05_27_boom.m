close all
clear

sample_date = '2019-05-27' ;
input_swcs_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/shared_tracing/Finished_Neurons/%s', sample_date) ;
input_sample_folder_path = sprintf('/nrs/mouselight/SAMPLES/%s', sample_date) ;
input_fragments_file_path = ...
    sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/build-brain-output/frags.mat', sample_date) ;
intermediate_results_folder_path = ...
    sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/build-brain-output/extract-fragments-near-first-pass-intermediate-results', ...
            sample_date) ;
output_swcs_folder_path = ...
    sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/build-brain-output/updated-frags-within-45-um-of-first-pass-neurons', sample_date) ;

% Get the metadata for the rendered sample
[shape_xyz, origin, spacing, jaws_origin] = load_sample_shape_origin_and_spacing(input_sample_folder_path);

% Load the fragments from disk                                                
tic_id = tic() ;
fragments_as_swc_arrays_in_erhan_coords = load_anonymous(input_fragments_file_path) ;
elapsed_time = toc(tic_id) ;
fprintf('Elapsed time to load in all fragments from a single .mat was %g seconds\n', elapsed_time) ;

% Extract all the fragment nodes, along with a look-up table to map them back to
% the fragment they came from
fragment_count = length(fragments_as_swc_arrays_in_erhan_coords) ;
[xyz_in_erhan_coords_from_fragment_node_id, fragment_id_from_fragment_node_id] = point_cloud_from_forest(fragments_as_swc_arrays_in_erhan_coords) ;
fragment_node_count = length(fragment_id_from_fragment_node_id) ;

% Check that the fragment nodes are on the grid as we expect
%     all_fragment_ijks_unrounded = (all_fragment_xyzs_in_erhan_coords - origin) ./ spacing + 0.5 ;  % these should be one-based ijks
%     all_fragment_ijks = round(all_fragment_ijks_unrounded) ;  % these should be one-based ijks
%     fragment_offsets = all_fragment_ijks_unrounded - all_fragment_ijks ;
%     assert( all( all( abs(fragment_offsets) < spacing/1024 ) ) ) ;  % the offsets should just be floating-point error                              
fraction_of_uncorrected_fragment_nodes_at_voxel_centers = ...
    fraction_of_xyzs_at_voxel_centers_using_erhan_conventions(xyz_in_erhan_coords_from_fragment_node_id, spacing, origin)
assert( fraction_of_uncorrected_fragment_nodes_at_voxel_centers == 1) ;

% Convert fragment node xyz coords to voxel indices
unrounded_ijk1_from_fragment_node_id = (xyz_in_erhan_coords_from_fragment_node_id - origin) ./ spacing + 0.5 ;  % these should be one-based ijks
ijk1_from_fragment_node_id = round(unrounded_ijk1_from_fragment_node_id) ;  % these should be one-based ijks

% % Check the alignment of the fragments to the fluorescence signal
% desired_sample_count = 100 ;
% [f,a] = figure_to_check_alignment_to_stack(sample_folder_path, ijk1_from_fragment_node_id, desired_sample_count) ;

% Convert one-based voxel indices to voxel center coords in JaWS coord system
% all_uncorrected_fragment_xyzs = jaws_origin + spacing .* (all_uncorrected_fragment_ijks-1) ;  % um, n x 3
% all_fragment_xyzs = jaws_origin + spacing .* (all_fragment_ijks-1) ;  % um, n x 3
xyz_from_fragment_node_id = jaws_origin + spacing .* (ijk1_from_fragment_node_id-1) ;  % um, n x 3
    % ijk1_from_fragment_node_id are 1-based, we want to map them to voxel centers
    % in JaWS coordinates

% Load in the first-pass neurons                          
first_pass_neurons_and_names = compute_or_read_from_memo(intermediate_results_folder_path, ...
                                                         'first_pass_neurons_and_names',  ...
                                                         @()(collect_first_pass_neurons(input_swcs_folder_path))) ;
first_pass_neurons_as_swc_arrays =  first_pass_neurons_and_names.neurons ;
first_pass_neuron_names = first_pass_neurons_and_names.neuron_names ;

fragment_kd_tree  = KDTreeSearcher(xyz_from_fragment_node_id) ;        
blast_radius_in_um = 45 ;  % um
is_fragment_near_neurons_from_fragment_id = ...
    find_fragments_near_neurons(...
        first_pass_neurons_as_swc_arrays, ...
        first_pass_neuron_names, ...
        fragment_kd_tree, ...
        fragment_node_count, ...
        fragment_id_from_fragment_node_id, ...
        fragment_count, ...
        blast_radius_in_um) ;

% Filter the fragments, convert to jaws coords
nearby_fragments_in_erhan_coords = fragments_as_swc_arrays_in_erhan_coords(is_fragment_near_neurons_from_fragment_id) ;
nearby_fragments = convert_swc_arrays_to_jaws_coords(nearby_fragments_in_erhan_coords, origin, spacing, jaws_origin) ;
                                                                                    
%%
% Write the swc's that don't exist already                          
save_swc_arrays(output_swcs_folder_path, ...                         
                nearby_fragments) ;
