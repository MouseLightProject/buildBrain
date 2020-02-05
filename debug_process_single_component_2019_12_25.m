sample_date = '2019-12-25' ;
options.whole_brain_h5_p_map_file_path = ...
    sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/whole-brain-p-map-as-h5/whole-brain-p-map.h5', sample_date) ;
options.whole_brain_h5_p_map_dataset_path = '/prob0' ;
options.whole_brain_h5_p_map_properties_group_path = '/prob0_props' ;
% options.skelfolder = ...
%     sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/skeletonization', sample_date) ;

options.skeleton_graph_mat_file_path = ...
    sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/skeleton-graph.mat', sample_date) ;


options.output_folder_path = ...
    sprintf('/dev/null') ;

%options.writefull = 1 ;
%options.writefrag = 1 ;

options.do_visualize = false ;
%options.debug = 0 ;

options.size_threshold = 10 ;  % trees must have at least this many nodes in the undecimated version to be kept

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

%
% extract_full_trees_as_mats(options) ;
%

% Extract parameters from options
whole_brain_h5_p_map_file_path = options.whole_brain_h5_p_map_file_path ;
%whole_brain_h5_p_map_dataset_path = options.whole_brain_h5_p_map_dataset_path ;
whole_brain_h5_p_map_properties_group_path = options.whole_brain_h5_p_map_properties_group_path ;

%brainSize = h5parser(whole_brain_h5_p_map_file_path, whole_brain_h5_p_map_dataset_path) ;

origin_in_nm = h5read(whole_brain_h5_p_map_file_path, [whole_brain_h5_p_map_properties_group_path,'/origin']) ;
top_level_spacing_in_nm = h5read(whole_brain_h5_p_map_file_path, [whole_brain_h5_p_map_properties_group_path,'/spacing']) ;
levels_below_top_level = h5read(whole_brain_h5_p_map_file_path, [whole_brain_h5_p_map_properties_group_path,'/level']) ;
%params.outsiz = brainSize ;
%     ox = origin_in_nm(1) ;
%     oy = origin_in_nm(2) ;
%     oz = origin_in_nm(3) ;
sx = top_level_spacing_in_nm(1) ;
sy = top_level_spacing_in_nm(2) ;
sz = top_level_spacing_in_nm(3) ;
%params.level = levels_below_top_level ;

voxres = [sx sy sz]/2^(levels_below_top_level)/1e3 ;  % in um, at highest zoom level
%options.params = params ;   

% Break out the options structure
%params = options.params ;
skeleton_graph_mat_file_path = options.skeleton_graph_mat_file_path ;
output_folder_path = options.output_folder_path ;
size_threshold = options.size_threshold ;
length_threshold = options.length_threshold ;
do_visualize = options.do_visualize ;
maximum_core_count_desired = options.maximum_core_count_desired ;
do_force_computations = options.do_force_computations ;
do_all_computations_serially = options.do_all_computations_serially ;
%sizetheshold = options.sizethreshold ;
sampling_style = options.sampling_style ;
sampling_interval = options.sampling_interval ;
%length_threshold = options.length_threshold ;
largesampling = options.largesampling ;

% Break out the 'params'
%origin_in_nm = [ox oy oz] ;  % nm
%top_level_spacing_in_nm = [sx sy sz] ;
%levels_below_top_level = params.level ;
spacing_in_nm = top_level_spacing_in_nm / 2^levels_below_top_level ;
%voxres = params.voxres ;

% Create the output folder, if needed
%output_folder_path = options.output_folder_path ;
% if ~exist(output_folder_path, 'file') ,
%     mkdir(output_folder_path) ;
% end

% Load the skeleton graph
load(skeleton_graph_mat_file_path, 'skeleton_graph', 'skeleton_ijks') ;   

% "components" are sometimes called "connected components"
[component_from_raw_component_id, size_from_raw_component_id] = conncomp(skeleton_graph, 'OutputForm', 'cell') ;  
  % cell array, each element a 1d array of node ids in G
[component_from_component_id, size_from_component_id] = sort_components_by_size(component_from_raw_component_id, size_from_raw_component_id) ;

% fractionated_components_file_path = fullfile(output_folder_path, 'fractionated_components.mat') ;
% if exist(fractionated_components_file_path, 'file') && ~do_force_computations ,
%     load(fractionated_components_file_path, 'component_from_component_id', 'size_from_component_id', 'maximum_component_size') ;
% else
%     %maximum_component_size = inf ;
%     maximum_component_size = 10e6 ;
%     [component_from_component_id, size_from_component_id] = connected_components_with_fractionation(skeleton_graph, maximum_component_size) ;
%        % note that components are sort by size, largest first
%     save(fractionated_components_file_path, 'component_from_component_id', 'size_from_component_id', 'maximum_component_size', '-v7.3') ;
% end

component_count = length(component_from_component_id) ;
fprintf('Components are in hand!  There are %d of them.\n', component_count) ;
if component_count > 0 ,
    fprintf('The largest component contains %d nodes.\n', size_from_component_id(1)) ;
end

% This will be useful in various places, so do it once here    
A = skeleton_graph.adjacency ;  % adjacency matrix, node_count x node_count, with zeros on the diagonal

% Extract medium-sized components
is_medium_sized = (100<=size_from_component_id) & (size_from_component_id<=100) ;

% Don't process ones that are too small
%fprintf('Filtering out components smaller than %d nodes...\n', size_threshold) ;
component_id_from_processing_index = find(is_medium_sized) ;    
components_to_process_count = length(component_id_from_processing_index)  %#ok<NOPTS>
%fprintf('%d components left.\n', components_to_process_count) ;

component_id_from_will_process_index = component_id_from_processing_index ;

% Figure out how many components are 'big', we'll do those one at a
% time so as not to run out of memory   
%big_component_threshold = 1e4 ;
do_process_serially_from_will_process_index = true(size(component_id_from_will_process_index)) ;
component_id_from_serial_will_process_index = component_id_from_will_process_index(do_process_serially_from_will_process_index) ;
%component_id_from_parallel_will_process_index = component_id_from_will_process_index(~do_process_serially_from_will_process_index) ;
components_to_process_serially_count = length(component_id_from_serial_will_process_index) ;
%components_to_process_in_parallel_count = length(component_id_from_parallel_will_process_index) ;

%%
fprintf('Starting the serial for loop, going to process %d components...\n', components_to_process_serially_count) ;
%did_discard = false(component_count, 1) ;
parfor_progress(components_to_process_serially_count) ;
% Do the big ones in a regular for loop, since each requires a lot of
% memory    
%for process_serially_index = 2 ,
for process_serially_index = 1 : components_to_process_serially_count ,
    component_id = component_id_from_serial_will_process_index(process_serially_index) ;        
    % Process this component
    component = component_from_component_id{component_id} ;
    ijks_for_component = skeleton_ijks(component,:) ;        
    %G_for_component = G.subgraph(component) ;  % very very very slow!
    A_for_component = A(component, component) ;  
    process_single_component(output_folder_path, ...
                             component_id, ...
                             component, ...
                             component_count, ...
                             A_for_component, ...
                             ijks_for_component, ...
                             size_threshold, ...
                             length_threshold, ...
                             do_visualize, ...
                             origin_in_nm, ...
                             spacing_in_nm, ...
                             voxres, ...
                             sampling_style, ...
                             sampling_interval, ...
                             largesampling) ;

    % Update the progress bar
    parfor_progress() ;
end
parfor_progress(0) ;
