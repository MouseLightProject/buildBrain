function build_fragments_as_swcs_workflow1(configuration_file_path)
    options = configparser(configuration_file_path) ;
    
    input_folder_path = options.input_folder_path ;
    output_folder_path = options.output_folder_path ;
    maximum_core_count_desired = options.maximum_core_count_desired ;
    minimum_centerpoint_count_per_fragment = options.minimum_centerpoint_count_per_fragment ;
    
    % if ~isfield(options,'sampling')
    %     options.sampling = 'uni';
    % end
%     whole_brain_h5_p_map_file_path = options.whole_brain_h5_p_map_file_path ;
%     whole_brain_h5_p_map_dataset_path = options.whole_brain_h5_p_map_dataset_path ;
%     whole_brain_h5_p_map_properties_group_path = options.whole_brain_h5_p_map_properties_group_path ;
    
%     brain_size = h5parser(whole_brain_h5_p_map_file_path, whole_brain_h5_p_map_dataset_path) ;
% 
%     origin = h5read(whole_brain_h5_p_map_file_path, [whole_brain_h5_p_map_properties_group_path, '/origin']) ;
%     spacing = h5read(whole_brain_h5_p_map_file_path, [whole_brain_h5_p_map_properties_group_path, '/spacing']) ;
%     level = h5read(whole_brain_h5_p_map_file_path, [whole_brain_h5_p_map_properties_group_path, '/level']) ;
%     params.outsiz = brain_size ;
%     params.ox = origin(1) ;
%     params.oy = origin(2) ;
%     params.oz = origin(3) ;
%     params.sx = spacing(1) ;
%     params.sy = spacing(2) ;
%     params.sz = spacing(3) ;
%     params.level = level ;
% 
%     params.voxres = [params.sx params.sy params.sz]/2^(params.level)/1e3 ;  % in um
    %options.params = params;

    %[subs,~,A,~] = skel2graph(options) ;

    %Gin = graph(max(A,A')) ;
    workflow1_frags_as_swcs(input_folder_path, ...
                            output_folder_path, ...
                            maximum_core_count_desired, ...
                            minimum_centerpoint_count_per_fragment) ;
end
