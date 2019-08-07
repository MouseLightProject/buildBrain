function missing_edge_count = count_missing_edges(G, components)
    original_edge_count = size(G.Edges,1) 
    component_count = length(components) ;
    edges_per_component = zeros(1, component_count) ;
    for i = 1 : component_count ,
        component = components{i} ;
        G_component = G.subgraph(component) ;
        edges_per_component(i) = size(G_component.Edges,1) ;
    end
    edges_per_component
    missing_edge_count = original_edge_count - sum(edges_per_component)
end
