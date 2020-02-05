%  1
%  ^^
%  | \
%  2  3

child_parent_pairs = [ 2 1 ; 3 1]
dA = sparse(child_parent_pairs(:,1), child_parent_pairs(:,2), 1, 3, 3)
[branches_struct, branch_id_from_node_id] = getBranchesOld(dA)
branches_struct(1)  
branches_struct(2)
branches_struct(3)

% [branches_struct, branch_id_from_node_id] = getBranchesOld(dA)
% 
% branches_struct = 
% 
%   1×3 struct array with fields:
% 
%     set
%     parentnode
%     parentbranch
% 
% 
% branch_id_from_node_id =
% 
%      1     2     3
% 
% branches_struct(1)
% 
% ans = 
% 
%   struct with fields:
% 
%              set: [1×0 double]
%       parentnode: 1
%     parentbranch: 1
% 
% branches_struct(2)
% 
% ans = 
% 
%   struct with fields:
% 
%              set: 2
%       parentnode: 1
%     parentbranch: 1
% 
% branches_struct(3)
% 
% ans = 
% 
%   struct with fields:
% 
%              set: 3
%       parentnode: 1
%     parentbranch: 1
% 
% Note that the "root branch" contains zero nodes (judging by its .set field), its parent branch
% is itself, and its parent node is the root node.  But note that
% branch_id_from_node_id(1) == 1, suggesting the root node is in the root
% branch.  


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
[branches_struct, branch_id_from_node_id] = getBranchesOld(dA)
branches_struct(1)  
branches_struct(2)
branches_struct(3)
branches_struct(4)

% branches_struct = 
% 
%   1×4 struct array with fields:
% 
%     set
%     parentnode
%     parentbranch
% 
% 
% branch_id_from_node_id =
% 
%      1     2     3     3     4     4
% 
% 
% ans = 
% 
%   struct with fields:
% 
%              set: [1×0 double]
%       parentnode: 1
%     parentbranch: 1
% 
% 
% ans = 
% 
%   struct with fields:
% 
%              set: 2
%       parentnode: 1
%     parentbranch: 1
% 
% 
% ans = 
% 
%   struct with fields:
% 
%              set: [4 3]
%       parentnode: 2
%     parentbranch: 2
% 
% 
% ans = 
% 
%   struct with fields:
% 
%              set: [6 5]
%       parentnode: 2
%     parentbranch: 2
%
% Note that the "root branch" contains zero nodes, its parent branch
% is itself, and its parent node is the root node.  Note also that the root
% node is not in the .set field of any branch, but the
% branch_id_from_node_id result seems to indicate the root node is in the
% root branch.

%  1
%  ^^^
%  | \ \
%  2  3 4

child_parent_pairs = [ 2 1 ; ...
                       3 1 ; ...
                       4 1 ]
dA = sparse(child_parent_pairs(:,1), child_parent_pairs(:,2), 1, 4, 4)
[branches_struct, branch_id_from_node_id] = getBranchesOld(dA)
branches_struct(1)  
branches_struct(2)
branches_struct(3)
branches_struct(4)
branches_struct(5)

% branches_struct = 
% 
%   1×5 struct array with fields:
% 
%     set
%     parentnode
%     parentbranch
% 
% 
% branch_id_from_node_id =
% 
%      2     3     4     5
% 
% 
% ans = 
% 
%   struct with fields:
% 
%              set: [1×0 double]
%       parentnode: 1
%     parentbranch: 2
% 
% 
% ans = 
% 
%   struct with fields:
% 
%              set: [1×0 double]
%       parentnode: 1
%     parentbranch: 2
% 
% 
% ans = 
% 
%   struct with fields:
% 
%              set: 2
%       parentnode: 1
%     parentbranch: 2
% 
% 
% ans = 
% 
%   struct with fields:
% 
%              set: 3
%       parentnode: 1
%     parentbranch: 2
% 
% 
% ans = 
% 
%   struct with fields:
% 
%              set: 4
%       parentnode: 1
%     parentbranch: 2
% 

% This just seems like a bug.  The first two branches are identical.  And
% no node even claims to be in branch 1, according to branch_id_from_node_id.

