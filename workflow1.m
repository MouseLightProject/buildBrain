function workflow1(G,subs,opt)
    % Break out the options structure
    params = opt.params ;
    output_folder_path = opt.outfolder ;
    size_threshold = opt.sizethreshold ;
    length_threshold = opt.lengthThr ;
    do_visualize = opt.viz ;
    
    % Fill out opt if needed
    if ~isfield(opt,'medianFiltSize') ,
        opt.medianFiltSize = 10 ;
    end
    
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
    raw_color_from_component_id = jet(component_count);
    color_from_component_id = raw_color_from_component_id(randperm(component_count),:);
    SX = params.sx;
    SY = params.sy;
    SZ = params.sz;
    voxres = [SX SY SZ]/2^(params.level)/1e3; % in um
    params.voxres = voxres;
    runtic = tic;
    %Eout = [];
    %iter = 0;
    %%
    maximum_core_count_desired = 20 ;
    poolobj = gcp('nocreate');  % If no pool, do not create new one.
    if isempty(poolobj) ,
        parpool([1 maximum_core_count_desired]) ;
    end
    poolobj = gcp('nocreate');  % If no pool, do not create new one.
    core_count = poolobj.NumWorkers ;
    fprintf('Using %d cores.\n', core_count) ;
    do_skip = false(1,component_count) ;
    do_skip(size_from_component_id<=size_threshold) = true ;
    %%
    full_trees_output_folder_path = fullfile(output_folder_path, 'full') ;
    fragment_output_folder_path =  fullfile(output_folder_path, 'frags') ;
    if ~exist(full_trees_output_folder_path, 'dir') ,
        mkdir(full_trees_output_folder_path) ;
    end
    if ~exist(fragment_output_folder_path, 'dir') ,
        mkdir(fragment_output_folder_path) ;
    end
    extant_full_tree_file_names = simple_dir(full_trees_output_folder_path) ;
    is_initial_pass = isempty(extant_full_tree_file_names) ;
    %acc=0;
    if ~is_initial_pass ,
        parfor component_id = 1:component_count ,   %1:size(Y,2)%ib(2:500)%
            if size_from_component_id(component_id) > size_threshold ,  %& ismember(component_id,validC) %& ~any(component_id==skipThese)
                tree_swc_file_name = sprintf('auto_cc-%04d.swc',component_id) ;
                tree_swc_file_path = fullfile(full_trees_output_folder_path,tree_swc_file_name);
                if exist(tree_swc_file_path,'file')
                    do_skip(component_id) = true ;
                    %acc=acc+1;
                end
            else
                do_skip(component_id) = true ;
            end
        end
    end

    %% 
    %try parfor_progress(0);catch;end    
    component_ids_to_process = find(~do_skip) ;
    components_to_process = components(component_ids_to_process) ;
    components_to_process_count = length(component_ids_to_process) ;
    did_discard = false(components_to_process_count, 1) ;
    fprintf('Starting the big parfor loop, going to process %d components...\n', components_to_process_count) ;
    parfor_progress(components_to_process_count) ;
    parfor i = 1 : components_to_process_count ,
        % for each cluster run reconstruction
        %%
        component_id = component_ids_to_process(i) ;
        fprintf('Processing component with id %d...\n', component_id) ;
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
            fprintf('Component with id %d fails to meet the size threshold after pruning, so discarding.\n', component_id) ;
            did_discard(i) = true ;
            continue
        end
        %%
        inupdate = smoothtree(inupdate,opt);
        %%
        % [L,list] = getBranches(inupdate.dA);
    %     if 0
    %         outtree = inupdate;
    %     else
            % outtree_old = downSampleTree(inupdate,opt);
        outtree = sampleTree(inupdate, opt) ;
    %     end
        if do_visualize
            cla
            gplot3(inupdate.dA,[inupdate.X,inupdate.Y,inupdate.Z],'LineWidth',3);
            hold on
            gplot3(outtree.dA,[outtree.X,outtree.Y,outtree.Z],'--','LineWidth',3);
            drawnow
        end
        
        % WRITE full tree
        if opt.writefull ,
            write_full_tree_as_swc(full_trees_output_folder_path, component_id, outtree, params, color_from_component_id) ;
        end
        
        % WRITE fragments
        if opt.writefrag ,
            write_fragments_as_swcs(fragment_output_folder_path, component_id, outtree, opt) ;
        end
        
        % Update the progress bar
        parfor_progress() ;
    end
    parfor_progress(0) ;
    discarded_components_count = sum(did_discard) ;
    fprintf('Of %d processed components, %d were discarded.\n', components_to_process_count, discarded_components_count) ;
    toc(runtic) ;
end