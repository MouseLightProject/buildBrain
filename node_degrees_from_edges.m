function result = node_degrees_from_edges(edges, max_node_id)
    serial_edges = edges(:) ;
    if ~exist('max_node_id', 'var') || isempty(max_node_id) ,
        max_node_id = max(serial_edges) ;
    end    
    result = zeros(max_node_id, 1) ;
    for node_id = 1 : max_node_id ,
        result(node_id) = sum(serial_edges==node_id) ;        
    end
end
