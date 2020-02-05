load('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2019-12-25/skeleton-graph.mat', 'skeleton_graph', 'skeleton_ijks') ;

%%
[components, size_per_component] = conncomp(skeleton_graph,'OutputForm','cell') ;  % cell array, each element a 1d array of node ids in G
[sorted_components, sorted_size_per_component] = sort_components_by_size(components, size_per_component) ;
