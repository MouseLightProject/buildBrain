input_folder_path = ...
    '/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2019-10-04/build-brain-output/full-as-named-tree-mats-take-2' ;
output_folder_path = ...
    '/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2019-10-04/build-brain-output/frags-from-take-2-trees-with-2-plus-nodes-as-mats' ;       
maximum_core_count_desired = inf ;
minimum_centerpoint_count_per_fragment = 2 ;    
bounding_box_low_corner_xyz = [-inf -inf -inf] ;
bounding_box_high_corner_xyz = [+inf +inf +inf] ;

workflow1_frags_as_mats(input_folder_path, ...
                        output_folder_path, ...
                        maximum_core_count_desired, ...
                        minimum_centerpoint_count_per_fragment, ...
                        bounding_box_low_corner_xyz, ...
                        bounding_box_high_corner_xyz) ;
