function write_fragments_as_mat(fragment_output_folder_path, ...
                                component_id, ...
                                full_tree, ...
                                minimum_centerpoint_count_per_fragment, ...
                                bounding_box_low_corner_xyz, ...
                                bounding_box_high_corner_xyz)
    XYZ = [full_tree.X,full_tree.Y,full_tree.Z] ;  % um
    R = full_tree.R ;
    D = full_tree.D ;
    [L,~] = getBranches(full_tree.dA) ;
    maximum_fragment_count = length(L) ;
    
    % Generate colormap
    %raw_color_from_fragment_id = jet(maximum_fragment_count);
    %color_from_fragment_id = raw_color_from_fragment_id(randperm(maximum_fragment_count),:) ;
    
    % Generate the file name for this file
    fragments_mat_file_name = sprintf('auto-cc-%04d-fragments.mat', component_id) ;
    fragments_mat_file_path = fullfile(fragment_output_folder_path, fragments_mat_file_name);
    
    % If the output file already exists, skip it
    if exist(fragments_mat_file_path, 'file') ,
        return
    end
        
    % Pre-allocate array to hold the fragments
    fragments_as_swc_arrays = cell(1, maximum_fragment_count) ;
    
    % Generate all the fragments in-memory
    fragment_index = 0 ;
    for fragment_id = 1:maximum_fragment_count ,
        L_this_fragment = L(fragment_id) ;
        out = L_this_fragment.set ;
        if isempty(out)
            continue
        end                

        % Get the centerpoints for this fragment
        XYZ_this_fragment = XYZ(out,:);
        R_this_fragment = R(out);
        D_this_fragment = D(out);

        % Don't write an swc for fragments with too few centerpoints
        centerpoint_count = size(XYZ_this_fragment,1) ;
        if centerpoint_count < minimum_centerpoint_count_per_fragment ,
            continue
        end
        
        % % Center centerpoint coordinates on the center of mass
        fragment_centroid_xyz = mean(XYZ_this_fragment, 1) ;        
        % XYZ_this_fragment_centered = bsxfun(@minus, XYZ_this_fragment, fragment_centroid_xyz) ;
        
        % Don't write fragments with centroids outside the bounding box
        if all(bounding_box_low_corner_xyz <= fragment_centroid_xyz) && all(fragment_centroid_xyz <= bounding_box_high_corner_xyz) ,
            % do nothing, fragment centroid is within bouding box
        else
            continue
        end

        % Create .swc data matrix
        centerpoint_ids = (1:centerpoint_count)' ;
        parent_ids_except_first = (1:(centerpoint_count-1))' ;  % parent of centerpoint 2 is 1, parent of centerpoint 3 is 2, etc.
        parent_ids = vertcat(-1, parent_ids_except_first) ;  % .swc says the root should have parent == -1        
        fragment_as_swc_array = horzcat(centerpoint_ids, D_this_fragment, XYZ_this_fragment, R_this_fragment, parent_ids) ;
        
        % Bump the fragment index
        fragment_index = fragment_index + 1 ;
        
        % Add to cell array
        fragments_as_swc_arrays{fragment_index} = fragment_as_swc_array ;        
    end

    % trim the array
    fragment_count = fragment_index ;    
    fragments_as_swc_arrays = fragments_as_swc_arrays(1:fragment_count) ;

    % Save these fragments to a file
    save(fragments_mat_file_path, '-v7.3', 'fragments_as_swc_arrays') ;    
end
