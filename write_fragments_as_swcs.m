function write_fragments_as_swcs(fragment_output_folder_path, component_id, full_tree, minimum_centerpoint_count_per_fragment)
    XYZ = [full_tree.X,full_tree.Y,full_tree.Z] ;  % um
    R = full_tree.R ;
    D = full_tree.D ;
    [L,~] = getBranches(full_tree.dA) ;
    fragment_count = length(L) ;
    
    % Generate colormap
    raw_color_from_fragment_id = jet(fragment_count);
    color_from_fragment_id = raw_color_from_fragment_id(randperm(fragment_count),:) ;
    
    % Write each fragment to disk as a .swc file
    for fragment_id = 1:fragment_count ,
        L_this_fragment = L(fragment_id) ;
        out = L_this_fragment.set ;
        if isempty(out)
            continue
        end        
        
        % .swc file name
        fragment_swc_file_name = sprintf('auto-cc-%04d-branch-%04d.swc', component_id, fragment_id) ;
        fragment_swc_file_path = fullfile(fragment_output_folder_path, fragment_swc_file_name);

        % If the output file already exists, skip it
        if exist(fragment_swc_file_path, 'file') ,
            continue ;
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
        
        % Center centerpoint coordinates on the center of mass
        offset = mean(XYZ_this_fragment, 1) ;        
        XYZ_this_fragment_centered = bsxfun(@minus, XYZ_this_fragment, offset) ;

        % Create .swc data matrix
        centerpoint_ids = (1:centerpoint_count)' ;
        parent_ids_except_first = (1:(centerpoint_count-1))' ;  % parent of centerpoint 2 is 1, parent of centerpoint 3 is 2, etc.
        parent_ids = vertcat(-1, parent_ids_except_first) ;  % .swc says the root should have parent == -1        
        swcData = horzcat(centerpoint_ids, D_this_fragment, XYZ_this_fragment_centered, R_this_fragment, parent_ids) ;
        
        % Write swc
        fid = fopen(fragment_swc_file_path,'w');
        mytxt = sprintf('# Generated by pw skel algorithm\n');
        fprintf(fid,'%s',mytxt);
        mytxt = sprintf('# OFFSET %.6f %.6f %.6f\n',offset);
        fprintf(fid,'%s',mytxt);
        mytxt = sprintf('# COLOR %f,%f,%f\n', ...
                        color_from_fragment_id(fragment_id,1), ...
                        color_from_fragment_id(fragment_id,2), ...
                        color_from_fragment_id(fragment_id,3)) ;
        fprintf(fid,'%s',mytxt);
        mytxt = sprintf('# NAME %s\n',fragment_swc_file_name);
        fprintf(fid,'%s',mytxt);
        fprintf(fid,'%d %d %f %f %f %d %d\n',swcData');
        fclose(fid);
    end
end
