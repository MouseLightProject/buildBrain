trees_as_mat_folder_path = '/nrs/mouselight/cluster/classifierOutputs/2019-05-27/build-brain-output/full-as-mats' ;
trees_as_swc_folder_path = '/nrs/mouselight/cluster/classifierOutputs/2019-05-27/build-brain-output/full' ;

use_this_many_cores(48) ;
convert_full_trees_to_swc(trees_as_swc_folder_path, trees_as_mat_folder_path) ;
