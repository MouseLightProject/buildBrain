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

    % The main loop
    fprintf('Starting the big parfor loop, going to process %d components...\n', component_count) ;
    did_discard = false(component_count, 1) ;
    parfor_progress(components_to_process_count) ;
    for component_id = 1 : component_count ,
        % Extract this component
        component = component_from_component_id{component_id} ;
        component_size = size_from_component_id(component_id) ;
        
        % Output some info
        fprintf('Processing component with id %d, size is %d nodes...\n', component_id, component_size) ;
        
        % Check if too small, skip iter if so
        if component_size<=size_threshold ,
            % Update the progress bar
            did_discard(component_id) = true ;
            parfor_progress() ;
            continue
        end                   
        
        % Check if the file already exists, skip iter if so
        tree_mat_file_name = sprintf('auto-cc-%06d.mat', component_id) ;
        tree_mat_file_path = fullfile(output_folder_path, tree_mat_file_name);
        if exist(tree_mat_file_path, 'file') ,
            % Update the progress bar
            parfor_progress() ;
            continue
        end
        
        % Get the graph for this component
        subs_for_component = subs(component,:) ;
        G_for_component = G.subgraph(component) ;
        A_for_component = G_for_component.adjacency ;

        % Do something
        a_leaf_node_id = find(sum(A_for_component)==1, 1) ;
        [eout] = graphfuncs.buildgraph(A_for_component, a_leaf_node_id) ;
        inupdate.dA = sparse(eout(:,1),eout(:,2),1,component_size,component_size);
        inupdate.D = ones(component_size,1);
        inupdate.R = ones(component_size,1);
        inupdate.X = subs_for_component(:,1);
        inupdate.Y = subs_for_component(:,2);
        inupdate.Z = subs_for_component(:,3);
        %%
        deleteThese = NaN;
        while ~isempty(deleteThese) ,
            [inupdate, deleteThese] = prunTree(inupdate, length_threshold, voxres) ;
            if do_visualize
                hold on
                gplot3(inupdate.dA,[inupdate.X,inupdate.Y,inupdate.Z]);
                drawnow
            end
        end
        %%
        % shuffle root to one of the leafs for efficiency and not
        % splitting long stretches into half
        if size(inupdate.dA,1)>1
            [eoutprun] = graphfuncs.buildgraph(inupdate.dA);
            component_size = max(eoutprun(:));
            inupdate.dA = sparse(eoutprun(:,1),eoutprun(:,2),1,component_size,component_size);
        else
        end
        if do_visualize
            hold on
            gplot3(inupdate.dA,[inupdate.X,inupdate.Y,inupdate.Z],'LineWidth',3);
            drawnow
        end
        %%
        if length(inupdate.dA)<size_threshold
            %fprintf('Component with id %d fails to meet the size threshold after pruning, so discarding.\n', component_id) ;
            did_discard(component_id) = true ;
            continue
        end
        %%
        inupdate_smoothed = smoothtree(inupdate,options);
        %%
        % [L,list] = getBranches(inupdate_smoothed.dA);
    %     if 0
    %         outtree = inupdate_smoothed;
    %     else
            % outtree_old = downSampleTree(inupdate_smoothed,opt);
        outtree_in_voxels = sampleTree(inupdate_smoothed, options) ;
        outtree_in_voxels.units = 'voxels' ;
        if do_visualize
            cla
            gplot3(inupdate_smoothed.dA,[inupdate_smoothed.X,inupdate_smoothed.Y,inupdate_smoothed.Z],'LineWidth',3);
            hold on
            gplot3(outtree_in_voxels.dA,[outtree_in_voxels.X,outtree_in_voxels.Y,outtree_in_voxels.Z],'--','LineWidth',3);
            drawnow
        end        
        
        % Convert centerpoint from voxel coords to um
        outtree = convert_centerpoint_units_to_um(outtree_in_voxels, origin_in_nm, spacing_in_nm) ;
                
        % Write full tree as a .mat file
        save_tree_as_mat(tree_mat_file_path, component_id, outtree) ;
        
        % Update the progress bar
        parfor_progress() ;
    end
    parfor_progress(0) ;
    discarded_components_count = sum(did_discard) ;
    fprintf('Of %d processed components, %d were discarded b/c they were too small after pruning.\n', ...
            component_count, ...
            discarded_components_count) ;
    if components_to_process_count == discarded_components_count ,
        fprintf('This means you don''t need to run this function again --- all trees have been generated.\n')
    end
    toc(runtic) ;
end