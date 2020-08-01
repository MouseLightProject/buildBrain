% This is not the code used for the tree generation for 2018-10-01.
% This is becuase, much later, I want to see if all the branches in the
% trees are binary.

% Set parameters
sample_date = '2018-10-01' ;
whole_brain_h5_p_map_file_path = ...
    sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/whole-brain-p-map-as-h5/whole-brain-p-map.h5', sample_date) ;
whole_brain_h5_p_map_properties_group_path = '/prob0_props' ;
skeleton_graph_mat_file_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/skeleton-graph.mat', sample_date) ;
output_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/build-brain-output/full-as-named-tree-mats', sample_date) ;

size_threshold = 100 ;  % node count
smoothing_filter_width = 10 ;  % nodes
length_threshold = 10 ;  % um, leaf chains shorter than this are pruned off
do_visualize = false ;
sampling_interval = 5 ;  % um

maximum_core_count_desired = 20 ;
do_force_computations = false ;
do_all_computations_serially = false ;

% Run
build_full_trees_as_named_tree_mats
