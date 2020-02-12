function [node_ids_in_chain, end_node_id] = trace_chain_in_undirected(A, start_node_id)
    % Trace a chain from a start node to the next
    % nexus point (a nexus is a node with degree ~= 2).  
    %
    % On return, node_ids_in_chain contains the chain, in the order the
    % nodes were traversed.  end_node_id contains the node id that
    % terminated the chain, because either it had no edges besides the
    % just-traversed edge, or there was more than one edge besides the
    % just-traversed edge.  end_node_id is never a member of
    % node_ids_in_chain on return. If node_ids_in_chain is nonempty,
    % node_ids_in_chain(1)==start_node_id.
    %
    % Note that if the graph contains cycles, this may enter an
    % infinite loop.
    
    node_ids_in_chain = zeros(1,0) ;
    last_node_id_maybe = zeros(1,0) ;    
    current_node_id = start_node_id ;        
    is_done = false ;
    while ~is_done ,
        % Get the adjacent nodes
        adjacent_node_ids = find(A(current_node_id,:)) ;
        % Filter out the last node to get the set of possible next nodes
        next_node_ids = setdiff(adjacent_node_ids, last_node_id_maybe) ;
        % See if the set of next nodes is a singleton set
        if isscalar(next_node_ids) ,
            % Add the current node to the chain
            node_ids_in_chain = [node_ids_in_chain current_node_id] ; %#ok<AGROW>  % add the current node to the list            
            % Set of next nodes is singleton, so
            % set up for next iteration
            last_node_id_maybe = current_node_id ;
            current_node_id = next_node_ids ;
        else                    
            % If the set of next nodes is not a singleton set, then exit the
            % loop
            end_node_id = current_node_id ;
            is_done = true ;
        end
    end    
end
