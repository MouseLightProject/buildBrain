input_file_path = ...
    '/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2019-10-04/build-brain-output/all-frags-from-take-2-trees-with-2-plus-nodes.mat' ;
output_folder_path = ...
    '/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2019-10-04/build-brain-output/frags-from-take-2-trees-with-2-plus-nodes-sans-melanin-3' ;
maximum_core_count_desired = inf ;

% Get the pool ready
use_this_many_cores(maximum_core_count_desired) ;

% Load the starting fragments
tic_id = tic() ;
load(input_file_path, 'fragments_as_named_trees') ;
elapsed_time = toc(tic_id) ;
fragment_count = length(fragments_as_named_trees) ;
fprintf('Elapsed time to load in all %d fragments from single .mat was %g seconds\n', fragment_count, elapsed_time) ;

% Get the centroid of each fragment
centroid_xyz_from_fragment_index = zeros(fragment_count, 3) ;
for fragment_index = 1 : fragment_count ,    
    tree = fragments_as_named_trees(fragment_index) ;
    centroid_xyz_from_fragment_index(fragment_index,:) = mean(tree.xyz,1) ;
end

% Calculate the number of neighbors around each centroid
kd_tree = KDTreeSearcher(centroid_xyz_from_fragment_index) ;
radius = 30 ;  % um
fragment_count_within_radius_from_fragment_index = zeros(fragment_count, 1) ;
parfor_progress(fragment_count) ;
parfor fragment_index = 1 : fragment_count ,    
    test_point = centroid_xyz_from_fragment_index(fragment_index,:) ;
    neighbors_in_cell_array = kd_tree.rangesearch(test_point, radius) ;
    neighbor_fragment_indices = neighbors_in_cell_array{1} ;
    neighbor_count = length(neighbor_fragment_indices)-1 ;  % don't want to count self
    fragment_count_within_radius_from_fragment_index(fragment_index) = neighbor_count ;
    parfor_progress() ;
end
parfor_progress(0) ;

[counts_from_bin_index, edges_from_bin_index] = histcounts(fragment_count_within_radius_from_fragment_index) ;
bin_center_from_bin_index = (edges_from_bin_index(1:end-1) + edges_from_bin_index(2:end))/2

figure;
bar(bin_center_from_bin_index, counts_from_bin_index) ;

% Get the node count for each fragment
node_count_from_fragment_index = zeros(fragment_count, 1) ;
for fragment_index = 1 : fragment_count ,    
    tree = fragments_as_named_trees(fragment_index) ;
    node_count_from_fragment_index(fragment_index) = size(tree.xyz,1) ;
end

% Plot node count vs neighbor count
figure('color', 'w') ;
plot(node_count_from_fragment_index, fragment_count_within_radius_from_fragment_index, '.k')
xlabel('Node count') ;
ylabel('Neighbor count') ;

%
% -2:
% is_melanin_from_fragment_index = (node_count_from_fragment_index<=5) & (fragment_count_within_radius_from_fragment_index>10) ;
% -3:
is_melanin_from_fragment_index = (node_count_from_fragment_index<=25) & (fragment_count_within_radius_from_fragment_index>10) ;
centroid_xyz_from_melanin_fragment_index = centroid_xyz_from_fragment_index(is_melanin_from_fragment_index,:) ;
centroid_xyz_from_non_melanin_fragment_index = centroid_xyz_from_fragment_index(~is_melanin_from_fragment_index,:) ;

min_xyz = min(centroid_xyz_from_fragment_index, [], 1) ;
max_xyz = max(centroid_xyz_from_fragment_index, [], 1) ;


figure('color', 'w') ;
plot3(centroid_xyz_from_non_melanin_fragment_index(:,1), centroid_xyz_from_non_melanin_fragment_index(:,2), centroid_xyz_from_non_melanin_fragment_index(:,3), '.b')
hold on
plot3(centroid_xyz_from_melanin_fragment_index(:,1), centroid_xyz_from_melanin_fragment_index(:,2), centroid_xyz_from_melanin_fragment_index(:,3), '.r')
hold off
xlabel('x (um)') ;
ylabel('y (um)') ;
zlabel('z (um)') ;
xlim([min_xyz(1) max_xyz(1)]) ;
ylim([min_xyz(2) max_xyz(2)]) ;
zlim([min_xyz(3) max_xyz(3)]) ;
set(gca, 'DataAspectRatio', [1 1 1]) ;
set(gca, 'PlotBoxAspectRatioMode', 'manual') ;
set(gca, 'CameraViewAngleMode', 'manual') ;

figure('color', 'w') ;
plot(centroid_xyz_from_non_melanin_fragment_index(:,1), centroid_xyz_from_non_melanin_fragment_index(:,2), '.b')
hold on
plot(centroid_xyz_from_melanin_fragment_index(:,1), centroid_xyz_from_melanin_fragment_index(:,2), '.r')
hold off
xlabel('x (um)') ;
ylabel('y (um)') ;
xlim([min_xyz(1) max_xyz(1)]) ;
ylim([min_xyz(2) max_xyz(2)]) ;
set(gca, 'YDir', 'reverse') ;
set(gca, 'DataAspectRatio', [1 1 1]) ;

% Get the non-melanin fragments
non_melanin_fragments_as_named_trees = fragments_as_named_trees(~is_melanin_from_fragment_index) ;

% Output to a folder of .swcs
save_swcs_from_named_trees(output_folder_path, ...                         
                           non_melanin_fragments_as_named_trees)
     


