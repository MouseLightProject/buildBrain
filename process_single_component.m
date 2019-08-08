function process_single_component(output_folder_path, ...
                                  component_id, ...
                                  component, ...
                                  G, ...
                                  subs, ...
                                  size_threshold, ...
                                  length_threshold, ...
                                  do_visualize, ...
                                  origin_in_nm, ...
                                  spacing_in_nm, ...
                                  voxres, ...
                                  options)
    component_size= length(component) ;
    
    % Output some info
    fprintf('Processing component with id %d, size is %d nodes...\n', component_id, component_size) ;

    % Check if too small, return if so
    if component_size<=size_threshold ,
        return
    end                   

    % Check if the file already exists, return if so
    tree_mat_file_name = sprintf('auto-cc-%06d.mat', component_id) ;
    tree_mat_file_path = fullfile(output_folder_path, tree_mat_file_name);
    if exist(tree_mat_file_path, 'file') ,
        return
    end

    % Get the graph for this component
    subs_for_component = subs(component,:) ;
    G_for_component = G.subgraph(component) ;
    A_for_component = G_for_component.adjacency ;

    % Do something
    root_node_id = find(sum(A_for_component)==1, 1) ;
    [eout] = graphfuncs.buildgraph(A_for_component, root_node_id) ;
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
        inupdate_A = max(inupdate.dA, inupdate.dA') ;  % make into an undirected graph adjacency matrix
        eoutprun = graphfuncs.buildgraph(inupdate_A, root_node_id) ;
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
        return
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
end    
