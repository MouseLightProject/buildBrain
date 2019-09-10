function build_fragments_as_swcs_workflow1_2019_05_27_small_chunk()
    input_folder_path = '/nrs/mouselight/cluster/classifierOutputs/2019-05-27/build-brain-output/full-as-mats' ;
    output_folder_path = '/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2019-05-27/build-brain-output/frags-small-chunk' ;
    maximum_core_count_desired = inf ;
    minimum_centerpoint_count_per_fragment = 1 ;    
    bounding_box_low_corner_xyz  = [76577 -inf  33856] ;
    bounding_box_high_corner_xyz = [78321 16145 36261] ;
    
    workflow1_frags_as_swcs(input_folder_path, ...
                            output_folder_path, ...
                            maximum_core_count_desired, ...
                            minimum_centerpoint_count_per_fragment, ...
                            bounding_box_low_corner_xyz, ...
                            bounding_box_high_corner_xyz) ;                            
end
