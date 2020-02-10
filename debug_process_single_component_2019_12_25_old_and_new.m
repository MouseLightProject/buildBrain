sample_date = '2019-12-25' ;
whole_brain_h5_p_map_file_path = ...
    sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/whole-brain-p-map-as-h5/whole-brain-p-map.h5', sample_date) ;
whole_brain_h5_p_map_dataset_path = '/prob0' ;
whole_brain_h5_p_map_properties_group_path = '/prob0_props' ;
% options.skelfolder = ...
%     sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/skeletonization', sample_date) ;

skeleton_graph_mat_file_path = ...
    sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/skeleton-graph.mat', sample_date) ;


%options.writefull = 1 ;
%options.writefrag = 1 ;

options.viz = false ;
%options.debug = 0 ;

options.sizethreshold = 10 ;  % trees must have at least this many nodes in the undecimated version to be kept

% overlap of running windows on skeletonization
options.fullh = 15 ;

%    length_threshold = 10;% (um) lengthThr for each leaf branch to decide to prune
%    largesampling = 200;%(um) make it more sparse for larger segments
options.length_threshold = 10 ;
options.largesampling = 200 ;
% sampling options are 'uni' [default] for uniform sampling (for 3D) and 'curv' for curvature weighted sampling (for 2D)
options.sampling_style = 'uni' ;  % 'uni' or 'curv'
options.sampling_interval = 5 ;  % keep every this many nodes when downsampling trees

%post segmentation parameters
options.tag = '' ;
options.prune = 1 ;
options.filterGraph = 0 ;
options.graph2branch = 1 ;

options.maximum_core_count_desired = 40 ;

options.do_force_computations = false ;
options.do_all_computations_serially = false ;


brainSize = h5parser(whole_brain_h5_p_map_file_path, whole_brain_h5_p_map_dataset_path) ;

origin_in_nm = h5read(whole_brain_h5_p_map_file_path, [whole_brain_h5_p_map_properties_group_path,'/origin']) ;
top_level_spacing_in_nm = h5read(whole_brain_h5_p_map_file_path, [whole_brain_h5_p_map_properties_group_path,'/spacing']) ;
levels_below_top_level = h5read(whole_brain_h5_p_map_file_path, [whole_brain_h5_p_map_properties_group_path,'/level']) ;

params.outsiz = brainSize ;
params.ox = origin_in_nm(1) ;
params.oy = origin_in_nm(2) ;
params.oz = origin_in_nm(3) ;
params.sx = top_level_spacing_in_nm(1) ;
params.sy = top_level_spacing_in_nm(2) ;
params.sz = top_level_spacing_in_nm(3) ;
params.level = levels_below_top_level ;

params.voxres = [params.sx params.sy params.sz]/2^(params.level)/1e3 ; % in um



load(skeleton_graph_mat_file_path, 'skeleton_graph', 'skeleton_ijks') ;   

% "components" are sometimes called "connected components"
[component_from_raw_component_id, size_from_raw_component_id] = conncomp(skeleton_graph, 'OutputForm', 'cell') ;  
  % cell array, each element a 1d array of node ids in G
[component_from_component_id, size_from_component_id] = sort_components_by_size(component_from_raw_component_id, size_from_raw_component_id) ;



 
%%
forest_of_named_trees_old = named_trees_from_options_old(skeleton_graph, ...
                                                         skeleton_ijks, ...
                                                         component_from_component_id, ...
                                                         size_from_component_id, ...
                                                         options, ...
                                                         params)

%%
forest_of_named_trees = named_trees_from_options(skeleton_graph, ...
                                                 skeleton_ijks, ...
                                                 component_from_component_id, ...
                                                 size_from_component_id, ...
                                                 options, ...
                                                 params)
