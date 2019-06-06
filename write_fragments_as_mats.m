function write_fragments_as_mats(fragment_output_folder_path, component_id, outtree, options)
    XYZ = [outtree.X,outtree.Y,outtree.Z]-1 ; % -1 one is heuristic/bug
    R = outtree.R ;
    D = outtree.D ;
    [L,~] = getBranches(outtree.dA) ;
    fragment_count = length(L) ;
    
    % % Generate colormap
    % raw_color_from_fragment_id = jet(fragment_count);
    % color_from_fragment_id = raw_color_from_fragment_id(randperm(fragment_count),:) ;
    
    % Write each fragment to disk as a .swc file
    for fragment_id = 1:fragment_count ,
        out = L(fragment_id).set ;
        if isempty(out)
            continue
        end        
        XYZout = XYZ(out,:);
        Rout = R(out);
        Dout = D(out);

        % transform location
        XYZout_um = pix2um(options.params, XYZout) ; % center anisotropy to compansate imresize
        %nout = size(XYZout,1) ;
        %swcData = [(1:nout)' Dout XYZout Rout (0:nout-1)'] ;
        %swcData(1,7) = -1;
        %offset = mean(swcData(:,3:5), 1) ;
        
        % Write mat file
        fragment_mat_file_name = sprintf('auto-cc-%04d-branch-%04d.mat', component_id, fragment_id) ;
        fragment_mat_file_path = fullfile(fragment_output_folder_path, fragment_mat_file_name) ;
        save(fragment_mat_file_path, '-mat', '-v7.3', 'component_id', 'fragment_id', 'XYZout_um', 'Dout', 'Rout') ;        
    end
end
