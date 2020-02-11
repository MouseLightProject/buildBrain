function result = smooth_tree(tree_as_dA_struct, filter_width)
    %SMOOTHTREE Summary of this function goes here
    %
    % [OUTPUTARGS] = SMOOTHTREE(INPUTARGS) Explain usage here
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

    % $Author: base $	$Date: 2016/03/29 10:17:20 $	$Revision: 0.1 $
    % Copyright: HHMI 2016
    
    branches = get_branches(tree_as_dA_struct.dA);

    xyz_from_node_id = [tree_as_dA_struct.X tree_as_dA_struct.Y tree_as_dA_struct.Z];
    %R = result.R;
    %D = result.D;

    for branch_index = 1:length(branches) ,
        this_branch = branches(branch_index) ;        
        node_ids = this_branch.node_ids ;
        working_node_ids = node_ids(2:end-1) ;
        if isempty(working_node_ids)
            continue
        end
        z_for_this_branch = xyz_from_node_id(working_node_ids,3) ;
        new_z_for_this_branch = medfilt1(medfilt1(z_for_this_branch, filter_width), filter_width) ;        
        xyz_from_node_id(working_node_ids,3) = new_z_for_this_branch ;
    end

    % Sub in the modified z coords to the result
    result = tree_as_dA_struct ;
    result.X = xyz_from_node_id(:,1) ;
    result.Y = xyz_from_node_id(:,2) ;
    result.Z = xyz_from_node_id(:,3) ;
end
