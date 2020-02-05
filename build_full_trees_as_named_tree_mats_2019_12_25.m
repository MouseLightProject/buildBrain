sample_date = '2019-12-25' ;
options.whole_brain_h5_p_map_file_path = ...
    sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/whole-brain-p-map-as-h5/whole-brain-p-map.h5', sample_date) ;
options.whole_brain_h5_p_map_dataset_path = '/prob0' ;
options.whole_brain_h5_p_map_properties_group_path = '/prob0_props' ;
% options.skelfolder = ...
%     sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/skeletonization', sample_date) ;

options.skeleton_graph_mat_file_path = ...
    sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/skeleton-graph.mat', sample_date) ;


options.outfolder = ...
    sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/build-brain-output/full-as-named-tree-mats', sample_date) ;

options.writefull = 1 ;
options.writefrag = 1 ;

options.viz = 0 ;
options.debug = 0 ;

options.sizethreshold = 100 ;
% decide crrop size from data as function of chunksize%cropSize = 1000

% overlap of running windows on skeletonization
options.fullh = 15 ;

%    lengthThr = 10;% (um) lengthThr for each leaf branch to decide to prune
%    largesampling = 200;%(um) make it more sparse for larger segments
options.lengthThr = 10 ;
options.largesampling = 200 ;
% sampling options are 'uni' [default] for uniform sampling (for 3D) and 'curv' for curvature weighted sampling (for 2D)
options.sampling = 'uni' ;
options.samplingInterval = 5 ;

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
output_folder_path = options.outfolder ;
size_threshold = options.sizethreshold ;
length_threshold = options.lengthThr ;
do_visualize = options.viz ;
maximum_core_count_desired = options.maximum_core_count_desired ;
do_force_computations = options.do_force_computations ;
do_all_computations_serially = options.do_all_computations_serially ;

% Break out the 'params'
%origin_in_nm = [ox oy oz] ;  % nm
%top_level_spacing_in_nm = [sx sy sz] ;
%levels_below_top_level = params.level ;
spacing_in_nm = top_level_spacing_in_nm / 2^levels_below_top_level ;
%voxres = params.voxres ;

% Create the output folder, if needed
%output_folder_path = options.outfolder ;
if ~exist(output_folder_path, 'file') ,
    mkdir(output_folder_path) ;
end

% Load the skeleton graph
load(skeleton_graph_mat_file_path, 'skeleton_graph', 'skeleton_ijks') ;   

% "components" are sometimes called "connected components"
%node_count = height(G.Nodes) ;
%[components, size_from_component_id] = conncomp(G,'OutputForm','cell') ;  % cell array, each element a 1d array of node ids in G
fractionated_components_file_path = fullfile(output_folder_path, 'fractionated_components.mat') ;
if exist(fractionated_components_file_path, 'file') && ~do_force_computations ,
    load(fractionated_components_file_path, 'component_from_component_id', 'size_from_component_id', 'maximum_component_size') ;  %#ok<NASGU>
else
    %maximum_component_size = inf ;
    maximum_component_size = 10e6 ;
    [component_from_component_id, size_from_component_id] = connected_components_with_fractionation(skeleton_graph, maximum_component_size) ;
       % note that components are sort by size, largest first
    save(fractionated_components_file_path, 'component_from_component_id', 'size_from_component_id', 'maximum_component_size', '-v7.3') ;
end
component_count = length(component_from_component_id) ;
fprintf('Components are in hand!  There are %d of them.\n', component_count) ;
if component_count > 0 ,
    fprintf('The largest component contains %d nodes.\n', size_from_component_id(1)) ;
end

% This will be useful in various places, so do it once here    
A = skeleton_graph.adjacency ;  % adjacency matrix, node_count x node_count, with zeros on the diagonal

%A_lower = tril(A,-1) ;  % lower-triangular part of A
%S = length(CompsC);
% size_from_component_id = cellfun(@length, components) ;
%   % the "component id" of a component is the index of that component in
%   % components
%component_id_from_component_id = 1:component_count ;
%[ia,ib]=sort(Y,'descend');
% algorithmkeys = {'spb','dij','bel','spa'};
% algorithm = 2;
% debug_level = 0;
% directed = 0;
%W = [];
%raw_color_from_component_id = jet(component_count);
%color_from_component_id = raw_color_from_component_id(randperm(component_count),:);
%     SX = params.sx;
%     SY = params.sy;
%     SZ = params.sz;
%     voxres = [SX SY SZ]/2^(params.level)/1e3; % in um
%     params.voxres = voxres;
runtic = tic;
%Eout = [];
%iter = 0;
%%

% Create the output folder if it doesn't exist
if ~exist(output_folder_path, 'dir') ,
    mkdir(output_folder_path) ;
end

% Don't process ones that are too small
fprintf('Filtering out components smaller than %d nodes...\n', size_threshold) ;
is_too_small = (size_from_component_id<=size_threshold) ;
component_id_from_processing_index = find(~is_too_small) ;    
components_to_process_count = length(component_id_from_processing_index) ;    
fprintf('%d components left.\n', components_to_process_count) ;

% Figure out which ones already exist
fprintf('Filtering out components for which output already exists...\n') ;
progress_bar = progress_bar_object(components_to_process_count) ;
does_output_exist_from_processing_index = false(1, components_to_process_count) ;
if do_force_computations ,
    progress_bar.update(components_to_process_count) ;
else
    digits_needed_for_index = floor(log10(component_count)) + 1 ;
    tree_name_template = sprintf('auto-cc-%%0%dd', digits_needed_for_index) ;  % e.g. 'tree-%04d'      
    for processing_index = 1 : components_to_process_count ,
        component_id = ...
            component_id_from_processing_index(processing_index) ;
        tree_name = sprintf(tree_name_template, component_id) ;
        tree_mat_file_name = sprintf('%s.mat', tree_name) ;
        tree_mat_file_path = fullfile(output_folder_path, tree_mat_file_name);
        does_output_exist_from_processing_index(processing_index) = logical(exist(tree_mat_file_path, 'file')) ;
        progress_bar.update(processing_index) ;
    end    
end
component_id_from_will_process_index = component_id_from_processing_index(~does_output_exist_from_processing_index) ;
will_process_count = length(component_id_from_will_process_index) ;        
fprintf('%d components left.\n', will_process_count) ;

% Figure out how many components are 'big', we'll do those one at a
% time so as not to run out of memory   
big_component_threshold = 1e4 ;
if do_all_computations_serially ,
    do_process_serially_from_will_process_index = true(size(component_id_from_will_process_index)) ;
else
    component_size_from_will_process_index = size_from_component_id(component_id_from_will_process_index) ;
    do_process_serially_from_will_process_index = (component_size_from_will_process_index>=big_component_threshold) ;
end
component_id_from_serial_will_process_index = component_id_from_will_process_index(do_process_serially_from_will_process_index) ;
component_id_from_parallel_will_process_index = component_id_from_will_process_index(~do_process_serially_from_will_process_index) ;
%component_size_from_serial_will_process_index = component_size_from_component_id(component_id_from_serial_will_process_index) ;
%component_size_from_parallel_will_process_index = component_size_from_component_id(component_id_from_parallel_will_process_index) ;
components_to_process_serially_count = length(component_id_from_serial_will_process_index) ;
components_to_process_in_parallel_count = length(component_id_from_parallel_will_process_index) ;

fprintf('Starting the serial for loop, going to process %d components...\n', components_to_process_serially_count) ;
%did_discard = false(component_count, 1) ;
parfor_progress(components_to_process_serially_count) ;
% Do the big ones in a regular for loop, since each requires a lot of
% memory    
for process_serially_index = 1 : components_to_process_serially_count ,
%for component_id = component_id_from_serial_will_process_index ,
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
                             options) ;

    % Update the progress bar
    parfor_progress() ;
end
parfor_progress(0) ;

% Do the small ones in a parfor loop, since memory is less of an issue
% for them    
fprintf('Starting the parallel for loop, will process %d components...\n', components_to_process_in_parallel_count) ;
use_this_many_cores(maximum_core_count_desired) ;
parfor_progress(components_to_process_in_parallel_count) ;
component_from_component_id_as_parpool_constant = parallel.pool.Constant(component_from_component_id) ;
subs_as_parpool_constant = parallel.pool.Constant(skeleton_ijks) ;
A_as_parpool_constant = parallel.pool.Constant(A) ;
parfor process_in_parallel_index = 1 : components_to_process_in_parallel_count ,          
    component_id = component_id_from_parallel_will_process_index(process_in_parallel_index) ;
    % Get all the parpool constants
    component_from_component_id_local = component_from_component_id_as_parpool_constant.Value ;
    subs_local = subs_as_parpool_constant.Value ;
    A_local = A_as_parpool_constant.Value ;
    % Process this component
    component = component_from_component_id_local{component_id} ;
    ijks_for_component = subs_local(component,:) ;
    %G_for_component = G_local.subgraph(component) ;  % slow slow slow
    A_for_component = A_local(component, component) ;  
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
                             options) ;

    % Update the progress bar
    parfor_progress() ;
end    
parfor_progress(0) ;

%     discarded_components_count = sum(did_discard) ;
%     fprintf('Of %d processed components, %d were discarded b/c they were too small after pruning.\n', ...
%             component_count, ...
%             discarded_components_count) ;
%     if components_to_process_count == discarded_components_count ,
%         fprintf('This means you don''t need to run this function again --- all trees have been generated.\n')
%     end
toc(runtic) ;
