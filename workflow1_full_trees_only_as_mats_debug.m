function workflow1_full_trees_only_as_mats_debug(G, ijks, options)
    % Break out the options structure
    params = options.params ;
    output_folder_path = options.outfolder ;
    size_threshold = options.sizethreshold ;
    length_threshold = options.lengthThr ;
    do_visualize = options.viz ;
    maximum_core_count_desired = options.maximum_core_count_desired ;

    % Break out the 'params'
    origin_in_nm = [params.ox params.oy params.oz] ;  % nm
    top_level_spacing = [params.sx params.sy params.sz] ;
    levels_below_top_level = params.level ;
    spacing_in_nm = top_level_spacing / 2^levels_below_top_level ;
    voxres = params.voxres ;
    
    % "components" are sometimes called "connected components"
    %node_count = height(G.Nodes) ;
    %[components, size_from_component_id] = conncomp(G,'OutputForm','cell') ;  % cell array, each element a 1d array of node ids in G
    fractionated_components_file_path = fullfile(output_folder_path, 'fractionated_components.mat') ;
    if exist(fractionated_components_file_path, 'file') && ~options.do_force_computations ,
        load(fractionated_components_file_path, 'component_from_component_id', 'size_from_component_id', 'maximum_component_size') ;  %#ok<NASGU>
    else
        %maximum_component_size = inf ;
        maximum_component_size = 10e6 ;
        [component_from_component_id, size_from_component_id] = connected_components_with_fractionation(G, maximum_component_size) ;
           % note that components are sort by size, largest first
        save(fractionated_components_file_path, 'component_from_component_id', 'size_from_component_id', 'maximum_component_size', '-v7.3') ;
    end
    component_count = length(component_from_component_id) ;
    fprintf('Components are in hand!  There are %d of them.\n', component_count) ;
    if component_count > 0 ,
        fprintf('The largest component contains %d nodes.\n', size_from_component_id(1)) ;
    end
    
    % Find components that impinge upon ROI
    padded_substack_heckbert_origin_xyz = [ 71178.0098371887          18796.0954300673          33772.2489377972 ] ;
    padded_substack_heckbert_far_corner_xyz = [ 71377.8774662781          19195.9629759657          33972.5159157184 ] ;
    origin_in_um = origin_in_nm/1000 ;
    spacing_in_um = spacing_in_nm/1000 ;
    xyzs = origin_in_um + spacing_in_um .* (ijks-1) ;
    is_in_roi_from_node_id = all((padded_substack_heckbert_origin_xyz<=xyzs) & (xyzs<=padded_substack_heckbert_far_corner_xyz), 2) ;
    
    % Create component_id_from_node_id array
    node_count = size(xyzs,1) ;
    component_id_from_node_id = zeros(node_count, 1) ;
    for component_id = 1 : component_count , 
        node_ids_in_this_component = component_from_component_id{component_id} ;
        component_id_from_node_id(node_ids_in_this_component) = component_id ;
    end
    
    component_ids_in_roi = unique(component_id_from_node_id(is_in_roi_from_node_id)) 
    
    
    % This will be useful in various places, so do it once here    
    A = G.adjacency ;  % adjacency matrix, node_count x node_count, with zeros on the diagonal
    
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

    % For now, assume we'll process all components
    % We'll winnow this down progressively...
    %component_ids_to_process = (1:component_count) ;
    component_ids_to_process = component_ids_in_roi ;
    components_to_process_count = length(component_ids_to_process) ;  
    fprintf('%d components at start.\n', components_to_process_count) ;
    
    % Don't process ones that are too small
    fprintf('Filtering out components smaller than %d nodes...\n', size_threshold) ;
    is_too_small_from_component_id = (size_from_component_id<=size_threshold) ;
    component_ids_that_are_too_small = find(is_too_small_from_component_id) ;    
    component_ids_to_process = setdiff(component_ids_to_process, component_ids_that_are_too_small) ;
    components_to_process_count = length(component_ids_to_process)
    fprintf('%d components left.\n', components_to_process_count) ;
    
%     % Figure out which ones already exist
%     fprintf('Filtering out components for which output already exists...\n') ;
%     progress_bar = progress_bar_object(components_to_process_count) ;
%     is_already_processed_from_component_id = false(1, component_count) ;
%     if options.do_force_computations ,
%         progress_bar.update(components_to_process_count) ;
%     else
%         digits_needed_for_index = floor(log10(component_count)) + 1 ;
%         tree_name_template = sprintf('auto-cc-%%0%dd', digits_needed_for_index) ;  % e.g. 'tree-%04d'      
%         for i = 1 : components_to_process_count ,
%             component_id = component_ids_to_process(i) ;
%             tree_name = sprintf(tree_name_template, component_id) ;
%             tree_mat_file_name = sprintf('%s.mat', tree_name) ;
%             tree_mat_file_path = fullfile(output_folder_path, tree_mat_file_name);
%             is_already_processed_from_component_id(component_id) = logical(exist(tree_mat_file_path, 'file')) ;
%             progress_bar.update(i) ;
%         end    
%     end
%     component_ids_already_processed = find(is_already_processed_from_component_id) ;
%     component_ids_to_process = setdiff(component_ids_to_process, component_ids_already_processed) ;
%     components_to_process_count = length(component_ids_to_process) ;        
%     fprintf('%d components left.\n', components_to_process_count) ;
    
    % Figure out how many components are 'big', we'll do those one at a
    % time so as not to run out of memory   
    big_component_threshold = 1e4 ;
    if options.do_all_computations_serially ,
        do_process_serially = true(size(component_ids_to_process)) ;
    else
        component_size_from_will_process_index = size_from_component_id(component_ids_to_process) ;
        do_process_serially = (component_size_from_will_process_index>=big_component_threshold) ;
    end
    component_ids_to_process_serially = component_ids_to_process(do_process_serially) ;
    component_ids_to_process_in_parallel = component_ids_to_process(~do_process_serially) ;
    %component_size_from_serial_will_process_index = component_size_from_component_id(component_id_from_serial_will_process_index) ;
    %component_size_from_parallel_will_process_index = component_size_from_component_id(component_id_from_parallel_will_process_index) ;
    components_to_process_serially_count = length(component_ids_to_process_serially) ;
    components_to_process_in_parallel_count = length(component_ids_to_process_in_parallel) ;

    fprintf('Starting the serial for loop, going to process %d components...\n', components_to_process_serially_count) ;
    %did_discard = false(component_count, 1) ;
    %parfor_progress(components_to_process_serially_count) ;
    % Do the big ones in a regular for loop, since each requires a lot of
    % memory    
    for i = 1 : components_to_process_serially_count ,
    %for component_id = component_id_from_serial_will_process_index ,
        component_id = component_ids_to_process_serially(i) ;        
        % Process this component
        component = component_from_component_id{component_id} ;
        ijks_for_component = ijks(component,:) ;        
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
        %parfor_progress() ;
    end
    %parfor_progress(0) ;

    % Do the small ones in a parfor loop, since memory is less of an issue
    % for them    
    fprintf('Starting the parallel for loop, will process %d components...\n', components_to_process_in_parallel_count) ;
    if components_to_process_in_parallel_count>0 ,
        use_this_many_cores(maximum_core_count_desired) ;
        parfor_progress(components_to_process_in_parallel_count) ;
        component_from_component_id_as_parpool_constant = parallel.pool.Constant(component_from_component_id) ;
        subs_as_parpool_constant = parallel.pool.Constant(ijks) ;
        A_as_parpool_constant = parallel.pool.Constant(A) ;
        parfor i = 1 : components_to_process_in_parallel_count ,          
            component_id = component_ids_to_process_in_parallel(i) ;
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
    end
    
%     discarded_components_count = sum(did_discard) ;
%     fprintf('Of %d processed components, %d were discarded b/c they were too small after pruning.\n', ...
%             component_count, ...
%             discarded_components_count) ;
%     if components_to_process_count == discarded_components_count ,
%         fprintf('This means you don''t need to run this function again --- all trees have been generated.\n')
%     end
    toc(runtic) ;
end
