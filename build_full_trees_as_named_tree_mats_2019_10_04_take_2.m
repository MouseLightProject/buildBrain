sample_date = '2019-10-04' ;
options.whole_brain_h5_p_map_file_path = ...
    sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/whole-brain-p-map-as-h5/whole-brain-p-map.h5', sample_date) ;
options.whole_brain_h5_p_map_dataset_path = '/prob0' ;
options.whole_brain_h5_p_map_properties_group_path = '/prob0_props' ;
options.skelfolder = ...
    sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/skeletonization', sample_date) ;

options.outfolder = ...
    sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/build-brain-output/full-as-named-tree-mats-take-2', sample_date) ;

options.writefull = 1 ;
options.writefrag = 1 ;

options.viz = 0 ;
options.debug = 0 ;

options.sizethreshold = 10 ;  % node count
% decide crrop size from data as function of chunksize%cropSize = 1000

% overlap of running windows on skeletonization
options.fullh = 15 ;

%    lengthThr = 10;% (um) lengthThr for each leaf branch to decide to prune
%    largesampling = 200;%(um) make it more sparse for larger segments
options.lengthThr = 10 ;  % um
options.largesampling = 200 ;
% sampling options are 'uni' [default] for uniform sampling (for 3D) and 'curv' for curvature weighted sampling (for 2D)
options.sampling = 'uni' ;
options.samplingInterval = 5 ;

%post segmentation parameters
options.tag = '' ;
options.prune = 1 ;
options.filterGraph = 0 ;
options.graph2branch = 1 ;

options.maximum_core_count_desired = 40 ;

options.do_force_computations = false ;
options.do_all_computations_serially = false ;

build_full_trees_as_mats_workflow1(options) ;
