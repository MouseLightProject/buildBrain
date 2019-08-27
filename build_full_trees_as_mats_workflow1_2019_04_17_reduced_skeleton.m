% make sure to add some header
sample_date = '2019-04-17' ;
this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;

options = struct() ;
options.whole_brain_h5_p_map_file_path = '/nrs/mouselight/cluster/classifierOutputs/2019-04-17/whole-brain-p-map-take-2/2019-04-17-whole-brain-p-map.h5' ;
options.whole_brain_h5_p_map_dataset_path = '/prob0' ;
options.whole_brain_h5_p_map_properties_group_path = '/prob0_props' ;
%options.skelfolder = '/nrs/mouselight/cluster/classifierOutputs/2019-04-17/skeletonization' ;
options.skelfolder = fullfile(this_folder_path, sprintf('%s-reduced-skeleton', sample_date)) ;

%options.outfolder = '/nrs/mouselight/cluster/classifierOutputs/2019-04-17/build-brain-output-staged-version' ;
options.outfolder = fullfile(this_folder_path, sprintf('%s-reduced-skeleton-output', sample_date)) ;

options.writefull = 1 ;
options.writefrag = 1 ;

options.viz = 0 ;
options.debug = 0 ;

options.sizethreshold = 100 ;
% decide crrop size from data as function of chunksize%cropSize = 1000

% overlap of running windows on skeletonization
options.fullh = 15 ;

%    lengthThr = 10;% (um) lengthThr for each leaf branch to decide to prune
%    largesampling = 200;%(um) make it more sparse for larger segments
options.lengthThr = 10 ;
options.largesampling = 200 ;
% sampling options are 'uni' [default] for uniform sampling (for 3D) and 'curv' for curvature weighted sampling (for 2D)
options.sampling = 'uni' ;
options.samplingInterval = 10 ;

%post segmentation parameters
options.tag = '' ;
options.prune = 1 ;
options.filterGraph = 0 ;
options.graph2branch = 1 ;

options.maximum_core_count_desired = 1 ;
options.do_force_computations = true ;
options.do_all_computations_serially = true ;

% Call the main function
build_full_trees_as_mats_workflow1(options) ;
