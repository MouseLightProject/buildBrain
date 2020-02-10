function result = named_trees_from_options(skeleton_graph, ...
                                           skeleton_ijks, ...
                                           component_from_component_id, ...
                                           size_from_component_id, ...
                                           options, ...
                                           params)
    % break out options, also populate the options.params field
%     whole_brain_h5_p_map_file_path = options.whole_brain_h5_p_map_file_path ;
%     whole_brain_h5_p_map_dataset_path = options.whole_brain_h5_p_map_dataset_path ;
%     whole_brain_h5_p_map_properties_group_path = options.whole_brain_h5_p_map_properties_group_path ;
% 
%     brainSize = h5parser(whole_brain_h5_p_map_file_path, whole_brain_h5_p_map_dataset_path) ;
% 
%     origin_in_nm = h5read(whole_brain_h5_p_map_file_path, [whole_brain_h5_p_map_properties_group_path,'/origin']) ;
%     top_level_spacing_in_nm = h5read(whole_brain_h5_p_map_file_path, [whole_brain_h5_p_map_properties_group_path,'/spacing']) ;
%     levels_below_top_level = h5read(whole_brain_h5_p_map_file_path, [whole_brain_h5_p_map_properties_group_path,'/level']) ;
%     params.outsiz = brainSize ;
%     params.ox = origin_in_nm(1) ;
%     params.oy = origin_in_nm(2) ;
%     params.oz = origin_in_nm(3) ;
%     params.sx = top_level_spacing_in_nm(1) ;
%     params.sy = top_level_spacing_in_nm(2) ;
%     params.sz = top_level_spacing_in_nm(3) ;
%     params.level = levels_below_top_level ;
% 
%     params.voxres = [params.sx params.sy params.sz]/2^(params.level)/1e3 ; % in um
%     options.params = params ;



    % % Extract parameters from options
    % whole_brain_h5_p_map_file_path = options.whole_brain_h5_p_map_file_path ;
    % %whole_brain_h5_p_map_dataset_path = options.whole_brain_h5_p_map_dataset_path ;
    % whole_brain_h5_p_map_properties_group_path = options.whole_brain_h5_p_map_properties_group_path ;
    % 
    % %brainSize = h5parser(whole_brain_h5_p_map_file_path, whole_brain_h5_p_map_dataset_path) ;
    % 
    % origin_in_nm = h5read(whole_brain_h5_p_map_file_path, [whole_brain_h5_p_map_properties_group_path,'/origin']) ;
    % top_level_spacing_in_nm = h5read(whole_brain_h5_p_map_file_path, [whole_brain_h5_p_map_properties_group_path,'/spacing']) ;
    % levels_below_top_level = h5read(whole_brain_h5_p_map_file_path, [whole_brain_h5_p_map_properties_group_path,'/level']) ;
    % %params.outsiz = brainSize ;
    % %     ox = origin_in_nm(1) ;
    % %     oy = origin_in_nm(2) ;
    % %     oz = origin_in_nm(3) ;
    % % sx = top_level_spacing_in_nm(1) ;
    % % sy = top_level_spacing_in_nm(2) ;
    % % sz = top_level_spacing_in_nm(3) ;
    % %params.level = levels_below_top_level ;
    % 
    top_level_spacing_in_nm = [params.sx params.sy params.sz] ;
    origin_in_nm = [params.ox params.oy params.oz] ;
    levels_below_top_level = params.level ;
    
    spacing_at_full_zoom = top_level_spacing_in_nm/2^(levels_below_top_level)/1e3 ;  % in um, at highest zoom level
    spacing_in_nm = top_level_spacing_in_nm / 2^levels_below_top_level ;
    % %options.params = params ;   
    % 
    % % Break out the options structure
    % %params = options.params ;
    %skeleton_graph_with_components_mat_file_path = options.skeleton_graph_with_components_mat_file_path ;
    % output_folder_path = options.output_folder_path ;
    size_threshold = options.sizethreshold ;
    length_threshold = options.length_threshold ;
    do_visualize = options.viz ;
    % maximum_core_count_desired = options.maximum_core_count_desired ;
    % do_force_computations = options.do_force_computations ;
    % do_all_computations_serially = options.do_all_computations_serially ;
    % %sizetheshold = options.sizethreshold ;
    sampling_style = options.sampling_style ;
    sampling_interval = options.sampling_interval ;
    %largesampling = options.largesampling ;

    % Break out the 'params'
    %origin_in_nm = [ox oy oz] ;  % nm
    %top_level_spacing_in_nm = [sx sy sz] ;
    %levels_below_top_level = params.level ;
    %voxres = params.voxres ;

    % Create the output folder, if needed
    %output_folder_path = options.output_folder_path ;
    % if ~exist(output_folder_path, 'file') ,
    %     mkdir(output_folder_path) ;
    % end

    % Load the skeleton graph
    %load(skeleton_graph_mat_file_path, 'skeleton_graph', 'skeleton_ijks') ;   
    %load(skeleton_graph_with_components_mat_file_path, 'skeleton_graph', 'skeleton_ijks', 'component_from_component_id', 'size_from_component_id') ;   

%     % "components" are sometimes called "connected components"
%     [component_from_raw_component_id, size_from_raw_component_id] = conncomp(skeleton_graph, 'OutputForm', 'cell') ;  
%       % cell array, each element a 1d array of node ids in G
%     [component_from_component_id, size_from_component_id] = sort_components_by_size(component_from_raw_component_id, size_from_raw_component_id) ;

    component_count = length(component_from_component_id) ;
    fprintf('There are %d components.\n', component_count) ;
    if component_count > 0 ,
        fprintf('The largest component contains %d nodes.\n', size_from_component_id(1)) ;
    end

    % This will be useful in various places, so do it once here    
    A = skeleton_graph.adjacency ;  % adjacency matrix, node_count x node_count, with zeros on the diagonal

    % Extract only medium-sized components
    is_medium_sized = (100<=size_from_component_id) & (size_from_component_id<=100) ;
    component_from_component_index_2 = component_from_component_id(is_medium_sized) ;
    component_id_from_component_id = (1:component_count) ;
    component_id_from_compenent_index_2 = component_id_from_component_id(is_medium_sized) ;
    %size_from_component_index_2 = size_from_component_id(is_medium_sized) ;
    component_index_2_count = length(component_from_component_index_2) ;
    
    %%
    max_iter_count = 50 ;
    iteration_count = min(max_iter_count, component_index_2_count) ;
    fprintf('Starting the for loop, going to process %d components...\n', iteration_count) ;
    parfor_progress(iteration_count) ;
    result = preallocate_forest_of_named_trees([iteration_count 1]) ;
    for component_index_2 = 1 : iteration_count ,
        % Process this component
        component_id = component_id_from_compenent_index_2(component_index_2) ;
        component = component_from_component_index_2{component_index_2} ;
        ijks_for_component = skeleton_ijks(component,:) ;        
        %G_for_component = G.subgraph(component) ;  % very very very slow!
        A_for_component = A(component, component) ;  
        named_tree = process_single_component_as_function(component_id, ...
                                                          component, ...
                                                          component_count, ...
                                                          A_for_component, ...
                                                          ijks_for_component, ...
                                                          size_threshold, ...
                                                          length_threshold, ...
                                                          do_visualize, ...
                                                          origin_in_nm, ...
                                                          spacing_in_nm, ...
                                                          spacing_at_full_zoom, ...
                                                          sampling_style, ...
                                                          sampling_interval) ;
        result(component_index_2) = named_tree ;

        % Update the progress bar
        parfor_progress() ;
    end
    parfor_progress(0) ;       
end
