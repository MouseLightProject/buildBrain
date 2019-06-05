this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;
configuration_file_path = fullfile(this_folder_path, 'config_files', 'config_buildBrain_20190417.cfg') ;
buildBrain_workflow1(configuration_file_path) ;
