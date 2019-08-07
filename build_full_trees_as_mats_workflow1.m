function build_full_trees_as_mats_workflow1(configuration_file_path)
    options = configparser(configuration_file_path);
    % if ~isfield(options,'sampling')
    %     options.sampling = 'uni';
    % end
    whole_brain_h5_p_map_file_path = options.whole_brain_h5_p_map_file_path ;
    whole_brain_h5_p_map_dataset_path = options.whole_brain_h5_p_map_dataset_path ;
    whole_brain_h5_p_map_properties_group_path = options.whole_brain_h5_p_map_properties_group_path ;
    
    brainSize = h5parser(whole_brain_h5_p_map_file_path, whole_brain_h5_p_map_dataset_path) ;

    origin = h5read(whole_brain_h5_p_map_file_path, [whole_brain_h5_p_map_properties_group_path,'/origin']) ;
    spacing = h5read(whole_brain_h5_p_map_file_path, [whole_brain_h5_p_map_properties_group_path,'/spacing']) ;
    level = h5read(whole_brain_h5_p_map_file_path, [whole_brain_h5_p_map_properties_group_path,'/level']) ;
    params.outsiz = brainSize ;
    params.ox = origin(1) ;
    params.oy = origin(2) ;
    params.oz = origin(3) ;
    params.sx = spacing(1) ;
    params.sy = spacing(2) ;
    params.sz = spacing(3) ;
    params.level = level ;
    
    params.voxres = [params.sx params.sy params.sz]/2^(params.level)/1e3 ; % in um
    options.params = params ;

    output_folder_path = options.outfolder ;
    graph_file_path = fullfile(output_folder_path, 'graph.mat') ;
    if exist(graph_file_path, 'file') ,
        load(graph_file_path, 'subs', 'Gin') ;
    else        
        [subs,~,A,~] = skel2graph(options) ;
        Gin = graph(max(A,A')) ;
        save(graph_file_path, 'subs', 'Gin', '-v7.3') ;
    end
    workflow1_full_trees_only_as_mats(Gin, subs, options) ;
end
