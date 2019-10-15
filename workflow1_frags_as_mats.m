function workflow1_frags_as_mats(input_folder_path, ...
                                 output_folder_path, ...
                                 maximum_core_count_desired, ...
                                 minimum_centerpoint_count_per_fragment, ...
                                 bounding_box_low_corner_xyz, ...
                                 bounding_box_high_corner_xyz)
    % Break out the options structure
    %params = options.params ;
    %output_folder_path = options.outfolder ;
    %size_threshold = options.sizethreshold ;
    %length_threshold = options.lengthThr ;
    %do_visualize = options.viz ;
    %maximum_core_count_desired =  opt.maximum_core_count_desired ;
    
    % "components" are sometimes called "connected components"
    %node_count = height(G.Nodes) ;
    %components = conncomp(G,'OutputForm','cell') ;  % cell array, each element a 1d array of node ids in G
    %component_count = length(components) ;
    %A = G.adjacency ;  % adjacency matrix, node_count x node_count, with zeros on the diagonal
    %A_lower = tril(A,-1) ;  % lower-triangular part of A
    %S = length(CompsC);
    %size_from_component_id = cellfun(@length, components) ;
      % the "component id" of a component is the index of that component in
      % components
    %component_id_from_component_id = 1:component_count ;
    %[ia,ib]=sort(Y,'descend');
    % algorithmkeys = {'spb','dij','bel','spa'};
    % algorithm = 2;
    % debug_level = 0;
    % directed = 0;
    %W = [];
    %raw_color_from_component_id = jet(component_count);
    %color_from_component_id = raw_color_from_component_id(randperm(component_count),:);
    %SX = params.sx;
    %SY = params.sy;
    %SZ = params.sz;
    %voxres = [SX SY SZ]/2^(params.level)/1e3; % in um
    %params.voxres = voxres;
    
    % Get the pool ready
    poolobj = gcp('nocreate');  % If no pool, do not create new one.
    if isempty(poolobj) ,
        parpool([1 maximum_core_count_desired]) ;
    end
    poolobj = gcp('nocreate');  % If no pool, do not create new one.
    core_count = poolobj.NumWorkers ;
    fprintf('Using %d cores.\n', core_count) ;
    
    %%
    %full_trees_folder_path = fullfile(output_folder_path, 'full-as-mats') ;
    %fragment_output_folder_path =  fullfile(output_folder_path, 'frags') ;
    if ~exist(output_folder_path, 'dir') ,
        mkdir(output_folder_path) ;
    end
    full_tree_file_names = simple_dir(fullfile(input_folder_path, 'auto-cc-*.mat')) ;
    full_tress_to_process_count = length(full_tree_file_names) ;
    tic_id = tic() ;
    fprintf('Starting the big parfor loop, going to process %d full trees...\n', full_tress_to_process_count) ;
    parfor_progress(full_tress_to_process_count) ;
    parfor full_tree_index = 1 : full_tress_to_process_count ,
        full_tree_file_name = full_tree_file_names{full_tree_index} ;
        full_tree_mat_file_path = fullfile(input_folder_path, full_tree_file_name) ;
        [component_id, outtree] = load_full_tree_from_mat(full_tree_mat_file_path) ;
        write_fragments_as_mat(output_folder_path, ...
                               component_id, ...
                               outtree, ...
                               minimum_centerpoint_count_per_fragment, ...
                               bounding_box_low_corner_xyz, ...
                               bounding_box_high_corner_xyz) ;
        
        % Update the progress bar
        parfor_progress() ;
    end
    parfor_progress(0) ;
    toc(tic_id) ;
end
