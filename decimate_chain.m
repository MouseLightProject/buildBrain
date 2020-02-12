function result = decimate_chain(node_ids_in_chain, xyz, desired_sampling_interval)
    % Get length of chain
    xyz_branch = xyz(node_ids_in_chain, :) ;
    dxyz_branch = diff(xyz_branch) ;        
    ds = sqrt( sum(dxyz_branch.^2, 2) ) ;  % length of each edge in the chain of nodes 
    s = [ 0 ; cumsum(ds) ] ;
    % Determine which nodes to keep, taking care that we always keep
    % the chain ends
    chain_length_in_um = s(end) ;
    if chain_length_in_um > desired_sampling_interval ,
        desired_s = (0:desired_sampling_interval:chain_length_in_um-desired_sampling_interval)' ;
        distance_matrix = pdist2(desired_s, s(:)) ;
        [~,idx] = min(distance_matrix, [], 2) ;
        result = [node_ids_in_chain(idx) node_ids_in_chain(end)] ;
    else
        result = [node_ids_in_chain(1) node_ids_in_chain(end)] ;
    end
end
