function build_fragments_as_swcs_workflow1_2019_05_27()
    this_file_path = mfilename('fullpath') ;
    this_folder_path = fileparts(this_file_path) ;
    configuration_file_path = fullfile(this_folder_path, 'config_files', 'build-fragments-2019-05-27.cfg') ;
    build_fragments_as_swcs_workflow1(configuration_file_path) ;
end
