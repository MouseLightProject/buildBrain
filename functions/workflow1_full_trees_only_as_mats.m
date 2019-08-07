function workflow1_full_trees_only_as_mats(G, subs, options)
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
    if exist(fractionated_components_file_path, 'file') ,
        load(fractionated_components_file_path, 'components', 'size_from_component_id', 'maximum_component_size') ;  %#ok<NASGU>
    else
        maximum_component_size = 10e6 ;
        [component_from_component_id, size_from_component_id] = connected_components_with_fractionation(G, maximum_component_size) ;
           % note that components are sort by size, largest first
        save(fractionated_components_file_path, 'components', 'size_from_component_id', 'maximum_component_size', '-v7.3') ;
    end
    component_count = length(component_from_component_id) ;
    %A = G.adjacency ;  % adjacency matrix, node_count x node_count, with zeros on the diagonal
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
    poolobj = gcp('nocreate');  % If no pool, do not create new one.
    if isempty(poolobj) ,
        parpool([1 maximum_core_count_desired]) ;
    end
    poolobj = gcp('nocreate');  % If no pool, do not create new one.
    core_count = poolobj.NumWorkers ;
    fprintf('Using %d cores.\n', core_count) ;
    
    % Create the output folder if it doesn't exist
    if ~exist(output_folder_path, 'dir') ,
        mkdir(output_folder_path) ;
    end

    % Figure out how many components are 'big', we'll do those one at a
    % time so as not to run out of memory   
    big_component_threshold = 1e6 ;
    big_component_count = sum(size_from_component_id>=big_component_threshold) ;
        
    % The main loop
    fprintf('Starting the big parfor loop, going to process %d components...\n', component_count) ;
    %did_discard = false(component_count, 1) ;
    parfor_progress(components_to_process_count) ;
    % Do the big ones in a regular for loop, since each requires a lot of
    % memory
    for component_id = 1 : big_component_count ,
        % Process this component
        component = component_from_component_id{component_id} ;
        process_single_component(component_id, component, G, subs, size_threshold, length_threshold, do_visualize, origin_in_nm, spacing_in_nm, voxres) ;
        
        % Update the progress bar
        parfor_progress() ;
    end

    % Do the small ones in a parfor loop, since memory is less of an issue
    % for them
    parfor component_id = big_component_count+1 : component_count ,
        % Process this component
        component = component_from_component_id{component_id} ;
        process_single_component(component_id, component, G, subs, size_threshold, length_threshold, do_visualize, origin_in_nm, spacing_in_nm, voxres) ;
        
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
end