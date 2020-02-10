function [output_dA_struct, did_prune_some] = prune_tree_old(input_dA_struct, length_threshold, spacing)
    % length_threshold should be in um.  Branches shorter than that will be
    % trimmed off.
    
    dA = input_dA_struct.dA;
    ijk1 = [input_dA_struct.X input_dA_struct.Y input_dA_struct.Z] ;  % node locations as voxel indices
    xyz = (ijk1-1) .* spacing ;  % these are relative to the origin, but we only use them for edge length computations
    %xyz = [input_dA_struct.X*spacing(1) input_dA_struct.Y*spacing(2) input_dA_struct.Z*spacing(3)];  % um
    branches = getBranches_old(dA);
    branch_count = length(branches);
    % find leaf branches
    node_indices_of_terminals = find(sum(dA)==0);
    %is_leaf_from_branch_index = zeros(1,branch_count);
    %deleteThese = [];
    node_count = size(xyz, 1) ;
    is_node_doomed = false(node_count, 1) ;
    for branch_index = 1:branch_count ,
        this_branch = branches(branch_index) ;
        node_indices_in_this_branch = this_branch.set;
        this_branch_parent_node_index = this_branch.parentnode ;
        if ~isempty(node_indices_in_this_branch) , ...
            if any(node_indices_of_terminals==node_indices_in_this_branch(1)) ,
                dr_to_parent = diff(xyz([node_indices_in_this_branch this_branch_parent_node_index],:)) ;
                length_to_parent = sqrt(sum(dr_to_parent.^2,2)) ;
                this_branch_length = sum(length_to_parent);
                if this_branch_length < length_threshold ,
                    %is_leaf_from_branch_index(branch_index) = 1;
                    %deleteThese = [deleteThese Liiset];
                    is_node_doomed(node_indices_in_this_branch) = true ;
                end
            end
        end
    end

    % Main output is the input minus the doomed nodes
    output_dA_struct = input_dA_struct;
    output_dA_struct.dA(is_node_doomed,:) = [];
    output_dA_struct.dA(:,is_node_doomed) = [];
    output_dA_struct.X(is_node_doomed) = [];
    output_dA_struct.Y(is_node_doomed) = [];
    output_dA_struct.Z(is_node_doomed) = [];
    output_dA_struct.R(is_node_doomed) = [];
    output_dA_struct.D(is_node_doomed) = [];

    % If requested, return whether we changed anything
    if nargout>=2 ,
        did_prune_some = any(is_node_doomed) ;
    end
end
