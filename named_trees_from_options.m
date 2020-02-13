function result = named_trees_from_options(skeleton_graph, ...
                                           skeleton_ijk1s, ...
                                           component_id_from_component_index, ...
                                           component_from_component_index, ...
                                           size_from_component_index, ...
                                           max_component_id, ...
                                           options, ...
                                           params)
    top_level_spacing_in_nm = [params.sx params.sy params.sz] ;
    origin_at_full_zoom = [params.ox params.oy params.oz]/1000 ;  % um
    levels_below_top_level = params.level ;
    
    spacing_at_full_zoom = top_level_spacing_in_nm/2^(levels_below_top_level)/1e3 ;  % in um, at highest zoom level
    size_threshold = options.sizethreshold ;
    smoothing_filter_width = size_threshold ; 
    length_threshold = options.length_threshold ;
    do_visualize = options.viz ;
    sampling_interval = options.sampling_interval ;

    component_count = length(component_from_component_index) ;
    fprintf('There are %d components.\n', component_count) ;
    if component_count > 0 ,
        fprintf('The largest component contains %d nodes.\n', size_from_component_index(1)) ;
    end

    % This will be useful in various places, so do it once here    
    A = skeleton_graph.adjacency ;  % adjacency matrix, node_count x node_count, with zeros on the diagonal
    skeleton_xyzs = origin_at_full_zoom + spacing_at_full_zoom .* (skeleton_ijk1s-0.5) ;  
      % NB: Using Erhan-style coords, here, to match the existing code.
      % Also, Erhan-style coords are Heckbertian, which I think are better
      % anyway.
    
    %%
    fprintf('Starting the for loop, going to process %d components...\n', component_count) ;
    parfor_progress(component_count) ;
    result = preallocate_forest_of_named_trees([component_count 1]) ;
    for component_index = 1 : component_count ,
        % Process this component
        component_id = component_id_from_component_index(component_index) ;
        component = component_from_component_index{component_index} ;
        xyzs_for_component = skeleton_xyzs(component,:) ;        
        %G_for_component = G.subgraph(component) ;  % very very very slow!
        A_for_component = A(component, component) ;  
        named_tree = process_single_component_as_function(component_id, ...
                                                          max_component_id, ...
                                                          A_for_component, ...
                                                          xyzs_for_component, ...
                                                          smoothing_filter_width, ...
                                                          length_threshold, ...
                                                          do_visualize, ...
                                                          sampling_interval) ;
        result(component_index) = named_tree ;

        % Update the progress bar
        parfor_progress() ;
    end
    parfor_progress(0) ;       
end
