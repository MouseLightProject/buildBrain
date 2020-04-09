sample_date = '2020-01-28-cameron-p-map' ;
trees_as_mat_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/build-brain-output/full-as-named-tree-mats', sample_date) ;
trees_as_swc_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/build-brain-output/full-as-swcs', sample_date) ;

use_this_many_cores(40) ;
convert_full_named_trees_as_mat_to_swc(trees_as_swc_folder_path, trees_as_mat_folder_path) ;
