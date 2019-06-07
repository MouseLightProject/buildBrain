function result = convert_centerpoint_units_to_um(tree, origin_in_nm, spacing_in_nm)
    if isequal(tree.units, 'um') ,
        result = tree ;
    elseif isequal(tree.units, 'voxels') ,
        result = tree ;
        r_in_voxels = [tree.X tree.Y tree.Z] ;
        r_in_um = um_from_voxels(r_in_voxels, origin_in_nm, spacing_in_nm) ;
        result.X = r_in_um(:,1) ;
        result.Y = r_in_um(:,2) ;
        result.Z = r_in_um(:,3) ;
        result.units = 'um' ;
    else
        error('Don''t know how to convert tree units to um') ;
    end
end
