function result = node_degrees_from_chains(chains, max_node_id)
    if ~exist('max_node_id', 'var') || isempty(max_node_id) ,
        max_node_id = max([chains{:}]) ;
    end
    result = zeros(max_node_id, 1) ;
    chain_count = length(chains) ;
    for i = 1 : chain_count , 
        chain = chains{i} ;
        head_id = chain(1) ;
        middle_ids = chain(2:end-1) ;
        tail_id = chain(end) ;
        result(head_id) = result(head_id) + 1 ;
        result(middle_ids) = result(middle_ids) + 2 ;
        result(tail_id) = result(tail_id) + 1 ;
    end
end
