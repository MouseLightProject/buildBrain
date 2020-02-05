function break_graph_into_chains(skeleton_graph, skeleton_ijk1s)
    [components, size_per_component] = conncomp(skeleton_graph,'OutputForm','cell') ;  % cell array, each element a 1d array of node ids in G
    [sorted_components, sorted_size_per_component] = sort_components_by_size(components, size_per_component) ;

end