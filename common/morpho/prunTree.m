function [output_dA_struct, did_prune_some] = prunTree(input_dA_struct, lengthThr, spacing)
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

    dA = input_dA_struct.dA;
    xyz = [input_dA_struct.X*spacing(1) input_dA_struct.Y*spacing(2) input_dA_struct.Z*spacing(3)];  % um
    L = getBranches(dA);
    numBranches = length(L);
    % find leaf branches
    termnodes = find(sum(dA)==0);
    leafbranches = zeros(1,numBranches);
    %deleteThese = [];
    node_count = size(dA, 1) ;
    is_node_doomed = false(node_count, 1) ;
    for ii=1:numBranches
        Liiset = L(ii).set;
        if ~isempty(Liiset) && any(termnodes==Liiset(1)) %& length(Liiset)<sizeThr
            lenBranch = sum(sqrt(sum(diff(xyz([Liiset L(ii).parentnode],:)).^2,2)));
            if lenBranch<lengthThr
                leafbranches(ii) = 1;
                %deleteThese = [deleteThese Liiset];
                is_node_doomed(Liiset) = true ;
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

    % If required, return a list of the deleted node_ids
    if nargout>=2 ,
        did_prune_some = any(is_node_doomed) ;
    end
end
