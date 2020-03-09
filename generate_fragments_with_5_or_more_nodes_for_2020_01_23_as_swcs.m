sample_date = '2020-01-23' ;
input_folder_path = ...
    sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/build-brain-output/full-as-named-tree-mats', ...
            sample_date) ;
output_folder_path = ...
    sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/build-brain-output/frags-with-5-or-more-nodes', ...
            sample_date) ;       
maximum_core_count_desired = inf ;
minimum_centerpoint_count_per_fragment = 5 ;    
bounding_box_low_corner_xyz = [-inf -inf -inf] ;
bounding_box_high_corner_xyz = [+inf +inf +inf] ;

generate_fragments_as_swcs(input_folder_path, ...
                           output_folder_path, ...
                           maximum_core_count_desired, ...
                           minimum_centerpoint_count_per_fragment, ...
                           bounding_box_low_corner_xyz, ...
                           bounding_box_high_corner_xyz) ;
