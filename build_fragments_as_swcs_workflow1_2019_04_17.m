this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;
configuration_file_path = fullfile(this_folder_path, 'config_files', 'config-build-brain-2019-04-17-take-2.cfg') ;
build_fragments_as_swcs_workflow1(configuration_file_path) ;
