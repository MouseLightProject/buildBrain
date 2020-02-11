function result = sample_tree(tree_as_dA_struct, spacing_at_full_zoom, desired_sampling_interval)
    %SAMPLETREE Summary of this function goes here
    %
    % [OUTPUTARGS] = SAMPLETREE(INPUTARGS) Explain usage here
    %
    % Inputs:
    %
    % Outputs:
    %
    % Examples:
    %
    % Provide sample usage code here
    %
    % See also: List related files here

    % $Author: base $	$Date: 2016/07/18 16:22:30 $	$Revision: 0.1 $
    % Copyright: HHMI 2016
    
    % Get a tree of chains, and the coords of each node
    branches = get_branches(tree_as_dA_struct.dA) ;
    ijk1 = [tree_as_dA_struct.X tree_as_dA_struct.Y tree_as_dA_struct.Z] ;
    
    %%
    % downsample branches
    branch_count = length(branches) ;
    downsampled_node_ids_from_branch_index = cell(branch_count, 1) ;
    for branch_index = 1 : branch_count ,
        %% 
        % Get a single branch
        branch = branches(branch_index) ;

        %% 
        % If this is the root branch, skip it
        if branch.parent_node_id == 0 ,
            continue
        end        

        %% 
        % Get the node ids for this branch
        branch_node_ids = branch.node_ids ;
        % branch_node_ids: +----x, includes distal node, excludes proximal node: indicies are
        % from distal to proximal
        % spacing is a function of branch length
        
        %%
        % Flip indicies so that lower indicies are more proximal.
        % Also, prepend the parent node.
        working_node_ids = [branch.parent_node_id flip(branch_node_ids)] ;
        
        % Get length of branch
        ijk1_branch = ijk1(working_node_ids,:) ;
        xyz_branch = ijk1_branch .* spacing_at_full_zoom ;  % Note that this is xyz with a funny origin, but that doesn't matter for now
        dxyz_branch = diff(xyz_branch) ;        
        ds = sqrt( sum(dxyz_branch.^2, 2) ) ;  % length of each edge in the chain of working nodes 
        s = [ 0 ; cumsum(ds) ] ;
        % Determine which nodes to keep, taking care that we always keep
        % the chain ends
        branch_length = s(end) ;
        if branch_length > desired_sampling_interval ,
            desired_s = (0:desired_sampling_interval:branch_length-desired_sampling_interval)' ;
            distance_matrix = pdist2(desired_s, s(:)) ;
            [~,idx] = min(distance_matrix, [], 2) ;
            downsampled_working_node_ids = [working_node_ids(idx) working_node_ids(end)] ;
        else
            downsampled_working_node_ids = [working_node_ids(1) working_node_ids(end)] ;
        end
        downsampled_node_ids_from_branch_index{branch_index} = downsampled_working_node_ids(2:end)' ;  % drop the proximal end
    end
    
    % Note that downsampled_node_ids_from_branch_index{branch_index} is
    % in proximal-to-distal order at this point.
    
    %
    % Now need to reassemble the downsampled branches into a dA struct
    %
    
    %%
    edges_from_branch_index = cell(branch_count, 1) ;
    for branch_index = 1 : branch_count ,
        inds = [ branches(branch_index).parent_node_id ; downsampled_node_ids_from_branch_index{branch_index}(:) ] ;
        to = inds(1:end-1);
        from = inds(2:end);
        edges_from_branch_index{branch_index} = [from(:) to(:)]';
    end
    raw_E = [edges_from_branch_index{:}]';

    %%
    [idx,~,ic] = unique(raw_E(:)) ;
    E = reshape(ic,[],2);
    dA = sparse(E(:,1),E(:,2),1,length(idx),length(idx));

    %%
    result.dA = dA;
    result.X = ijk1(idx,1);
    result.Y = ijk1(idx,2);
    result.Z = ijk1(idx,3);
    result.R = tree_as_dA_struct.R(idx) ;
    result.D = tree_as_dA_struct.D(idx) ;
end
