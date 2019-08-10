function [inupdate, did_prune_some] = prunTree(in_, lengthThr, res)
    %PRUNTREE Pruns a tree with a given length threshold
    % 
    % [OUTPUTARGS] = PRUNTREE(INPUTARGS) Explain usage here
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

    % $Author: base $	$Date: 2015/10/28 10:42:45 $	$Revision: 0.1 $
    % Copyright: HHMI 2015

    A = in_.dA;
    XYZ = [in_.X*res(1) in_.Y*res(2) in_.Z*res(3)];% in (um)
    L = getBranches(A);
    numBranches = length(L);
    % find leaf branches
    termnodes = find(sum(A)==0);
    leafbranches = zeros(1,numBranches);
    %deleteThese = [];
    node_count = size(A, 1) ;
    is_node_doomed = false(node_count, 1) ;
    for ii=1:numBranches
        Liiset = L(ii).set;
        if ~isempty(Liiset) && any(termnodes==Liiset(1)) %& length(Liiset)<sizeThr
            lenBranch = sum(sqrt(sum(diff(XYZ([Liiset L(ii).parentnode],:)).^2,2)));
            if lenBranch<lengthThr
                leafbranches(ii) = 1;
                %deleteThese = [deleteThese Liiset];
                is_node_doomed(Liiset) = true ;
            end
        end
    end

    % Main output is the input minus the doomed nodes
    inupdate = in_;
    inupdate.dA(is_node_doomed,:) = [];
    inupdate.dA(:,is_node_doomed) = [];
    inupdate.X(is_node_doomed) = [];
    inupdate.Y(is_node_doomed) = [];
    inupdate.Z(is_node_doomed) = [];
    inupdate.R(is_node_doomed) = [];
    inupdate.D(is_node_doomed) = [];

    % If required, return a list of the deleted node_ids
    if nargout>=2 ,
        did_prune_some = any(is_node_doomed) ;
    end
end
