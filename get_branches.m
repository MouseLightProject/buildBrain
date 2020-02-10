function [branches_struct, branch_id_from_node_id] = get_branches(dA)
    % Given directed adjacency matrix dA representing a rooted tree,
    % returns a tree of branches, where each "branch" is a chain from dA, and
    % which mimics the topology of the rooted tree dA.    
    %
    % Edges in dA should point rootward.
    % I.e. dA(i,j) == 1 implies that j is the next node rootward from i.    
    %
    % On return, branches_struct is branch_count x 1, with fields:
    %
    %     branches_struct(i).node_ids: the node ids in branch i, with the
    %         most-distal node first.  Each subsequent node id is the
    %         parent of the previous node it in the rooted tree represented
    %         by dA.
    %     branches_struct(i).parent_node_id: the node id of the parent node
    %         of branch i
    %     branches_struct(i).parent_branch_id: the branch id of the parent
    %         branch of branch i
    %
    % On return, branch_id_from_node_id is node_count x 1, and gives the
    % branch id of the branch containing each node.  Every node is a member
    % of some branch.
    %
    % Note that the root branch always contains only a single node, the
    % root node.  The parent_node_id of the root branch is 0 by convention,
    % and the parent_branch_id is 0 also.

%     % For debugging
%     G = digraph(dA) ;
%     figure() ;
%     plot(G) ;

    % Get the node count
    node_count = size(dA,1) ;
    
    % Find the root node
    out_degree_from_node_id = sum(dA,2) ;
    root_node_id = find(out_degree_from_node_id==0) ;  % the root node is the one with no outgoing edges
    
    % Find the branch nodes
    in_degree_from_node_id = sum(dA,1) ;
    branch_node_ids = find(in_degree_from_node_id>1) ;
    
    % Find the terminal node ids
    terminal_node_ids = find(in_degree_from_node_id==0) ;
    
    % from all branch and termnodes, get segmentsuntill another critical point
    % is found
    distal_tip_node_id_from_branch_index = unique([root_node_id branch_node_ids terminal_node_ids]) ;  
      % the root node might be a branch or terminal node, so do unique()
    branch_count = length(distal_tip_node_id_from_branch_index) ;
    
    % Make a look-up table for which nodes are branch tips
    is_branch_tip_from_node_id = false(node_count,1) ;
    is_branch_tip_from_node_id(distal_tip_node_id_from_branch_index) = true ;
    
    % Extract the parent look-up table from dA
    [~, parent_node_id_from_node_id] = graphtraverse(dA', root_node_id, 'DIRECTED', true) ;  % in this, the parent of the root is 0

    % Trace each branch from its distal end to its proximal end
    branches_struct = struct_with_shape_and_fields([branch_count 1], {'node_ids', 'parent_node_id', 'parent_branch_id'}) ;
    branch_id_from_node_id = zeros(1,node_count) ;
    for branch_id = 1 : branch_count ,
        %%
        distal_node_id = distal_tip_node_id_from_branch_index(branch_id) ;
        if distal_node_id == root_node_id ,
            branches_struct(branch_id).node_ids = root_node_id ;
            branches_struct(branch_id).parent_node_id = 0 ;
            branch_id_from_node_id(root_node_id) = branch_id ;
        else            
            [node_ids_in_this_branch, parent_of_this_branch_node_id] = ...
                trace_rootward(parent_node_id_from_node_id, is_branch_tip_from_node_id, distal_node_id) ;
            branches_struct(branch_id).node_ids = node_ids_in_this_branch ;
            branches_struct(branch_id).parent_node_id = parent_of_this_branch_node_id ;
            branch_id_from_node_id(node_ids_in_this_branch) = branch_id;
        end
    end
    
    % For each branch, determine the branch index of its parent branch
    for branch_id = 1 : branch_count ,
        parent_node_id = branches_struct(branch_id).parent_node_id ;
        if parent_node_id == 0 ,
            % this means this is the root branch
            branches_struct(branch_id).parent_branch_id = 0 ;
        else
            parent_branch_id = branch_id_from_node_id(parent_node_id) ;
            branches_struct(branch_id).parent_branch_id = parent_branch_id ;
        end
    end
end



function [node_ids_in_branch, parent_of_branch_node_id] = trace_rootward(parent_node_id_from_node_id, is_branch_tip_from_node_id, starting_node_id)
    % Trace a single branch rootward from a starting node to the next
    % branch point (or the root).  On return, node_ids_in_branch(1) is always
    % starting_node_id.  The node id of branch (or root) node that is
    % just-rootward of the branch is given by parent_of_branch_node_id.
    % Note that parent_of_branch_node_id is never an element of
    % node_ids_in_branch.
    %
    % On call, is_branch_tip_from_node_id should be true iff that node is a
    % 1) the root, 2) a branch node (has 3 or more neighbors in the undirected 
    % version of the graph), or 3) a terminal node.
    
    node_ids_in_branch = starting_node_id ;
    current_node_id = parent_node_id_from_node_id(starting_node_id) ;
    if current_node_id==0 ,
        % this means the starting node is the root, so can't trace further
    else
        while ~is_branch_tip_from_node_id(current_node_id) ,
            node_ids_in_branch = [node_ids_in_branch current_node_id] ; %#ok<AGROW>  % add the node to the list
            current_node_id = parent_node_id_from_node_id(current_node_id) ;
        end    
    end
    parent_of_branch_node_id = current_node_id ;
end
