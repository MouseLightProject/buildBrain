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
    components = conncomp(G,'OutputForm','cell') ;  % cell array, each element a 1d array of node ids in G
    component_count = length(components) ;
    A = G.adjacency ;  % adjacency matrix, node_count x node_count, with zeros on the diagonal
    A_lower = tril(A,-1) ;  % lower-triangular part of A
    %S = length(CompsC);
    size_from_component_id = cellfun(@length, components) ;
      % the "component id" of a component is the index of that component in
      % components
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

    % Figure out which components we can skip
    is_too_small = (size_from_component_id<=size_threshold) ;  % we'll skip these    
    is_already_done = false(1,component_count) ;        
    extant_full_tree_file_names = simple_dir(output_folder_path) ;    
    is_initial_pass = isempty(extant_full_tree_file_names) ;
    if ~is_initial_pass ,
        parfor component_id = 1:component_count ,
            if ~is_too_small(component_id) , 
                tree_mat_file_name = sprintf('auto-cc-%06d.mat', component_id) ;
                tree_mat_file_path = fullfile(output_folder_path, tree_mat_file_name);
                if exist(tree_mat_file_path, 'file') ,
                    is_already_done(component_id) = true ;
                end
            end
        end
    end
    do_skip = is_too_small | is_already_done ;
    do_process = ~do_skip ;

    %try parfor_progress(0);catch;end    
    component_ids_to_process = find(do_process) ;
    components_to_process = components(component_ids_to_process) ;
    components_to_process_count = length(component_ids_to_process) ;
    did_discard = false(components_to_process_count, 1) ;
    fprintf('Starting the big parfor loop, going to process %d components...\n', components_to_process_count) ;
    parfor_progress(components_to_process_count) ;
    parfor i = 1 : components_to_process_count ,
        % for each cluster run reconstruction
        %%
        component_id = component_ids_to_process(i) ;
        %fprintf('Processing component with id %d...\n', component_id) ;
        %component_id = component_id ;  % iter+1;
        component = components_to_process{i} ;  % find(Comps==component_id);
        subs_for_component = subs(component,:) ;  %#ok<PFBNS> % get back to matlab image coordinates
        component_size = length(component) ;
        % get lower portion to make it directed
        A_lower_for_component = A_lower(component,component) ;  %#ok<PFBNS> % faster
        %Gsub = G.subgraph(subidx);

        leafs = find(sum(A_lower_for_component,2)==0);%find(sum(max(Asub,Asub'))==1,1);

        [eout] = graphfuncs.buildgraph(A_lower_for_component,leafs(1));
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
            did_discard(i) = true ;
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
        
    %     end
                
        % Write full tree as a .mat file
        tree_mat_file_name = sprintf('auto-cc-%06d.mat', component_id) ;
        tree_mat_file_path = fullfile(output_folder_path, tree_mat_file_name);
        save_tree_as_mat(tree_mat_file_path, component_id, outtree) ;
        
        % Update the progress bar
        parfor_progress() ;
    end
    parfor_progress(0) ;
    discarded_components_count = sum(did_discard) ;
    fprintf('Of %d processed components, %d were discarded b/c they were too small after pruning.\n', ...
            components_to_process_count, ...
            discarded_components_count) ;
    if components_to_process_count == discarded_components_count ,
        fprintf('This means you don''t need to run this function again --- all trees have been generated.\n')
    end
    toc(runtic) ;
end