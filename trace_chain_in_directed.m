function [node_ids_in_chain, end_node_id] = trace_chain_in_directed(dA, start_node_id)
    % Trace a chain from a start node to the next
    % nexus point (a nexus is a node with degree ~= 2).  
    %
    % On return, node_ids_in_chain contains the chain, in the order
    % tranversed.  end_node_id contains the node id that terminated the
    % chain, because either it had no outgoing edges, or there was more
    % than one ongoing edge. end_node_id is ever a member of
    % node_ids_in_chain on return. If node_ids_in_chain is nonempty,
    % node_ids_in_chain(1)==start_node_id.
    %
    % Note that if the graph contains directed cycles, this may enter an
    % infinite loop.
    
    node_ids_in_chain = zeros(1,0) ;
    current_node_id = start_node_id ;
    is_done = false ;
    while ~is_done ,
        % Get the adjacent nodes
        adjacent_node_ids = find(dA(current_node_id,:)) ;
        % See if the set of next nodes is a singleton set
        if isscalar(adjacent_node_ids) ,
            % Add the current node to the chain
            node_ids_in_chain = [node_ids_in_chain current_node_id] ; %#ok<AGROW>  % add the current node to the list
            % Set of next nodes is singleton, so
            % set up for next iteration
            current_node_id = adjacent_node_ids ;
        else                    
            % If the set of next nodes is not a singleton set, then exit the
            % loop
            end_node_id = current_node_id ;
            is_done = true ;
        end
    end
end
