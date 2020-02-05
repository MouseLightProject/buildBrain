%  1
%  ^^
%  | \
%  2  3

child_parent_pairs = [ 2 1 ; 3 1]
dA = sparse(child_parent_pairs(:,1), child_parent_pairs(:,2), 1, 3, 3)
[branches_struct, branch_id_from_node_id] = get_branches(dA)
branches_struct(1)  
branches_struct(2)
branches_struct(3)



%  1
%  ^
%  | 
%  2  
%  ^^
%  | \
%  3  5
%  ^  ^
%  |  |
%  4  6

child_parent_pairs = [ 2 1 ; ...
                       3 2 ; ...
                       4 3 ; ...
                       5 2 ; ...
                       6 5 ] ;
dA = sparse(child_parent_pairs(:,1), child_parent_pairs(:,2), 1, 6, 6)
[branches_struct, branch_id_from_node_id] = get_branches(dA)
branches_struct(1)  
branches_struct(2)
branches_struct(3)
branches_struct(4)

%  1
%  ^^^
%  | \ \
%  2  3 4

child_parent_pairs = [ 2 1 ; ...
                       3 1 ; ...
                       4 1 ]
dA = sparse(child_parent_pairs(:,1), child_parent_pairs(:,2), 1, 4, 4)
[branches_struct, branch_id_from_node_id] = get_branches(dA)
branches_struct(1)  
branches_struct(2)
branches_struct(3)
branches_struct(4)

