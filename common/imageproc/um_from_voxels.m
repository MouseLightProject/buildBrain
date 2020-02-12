function result = um_from_voxels(r_in_voxels, origin_in_nm, spacing_in_nm)
    % This like pix2um, but more to ALT's liking.
    % r_in_voxels is n x 3, in x, y, z order, in units of voxel spacing (thus pure numbers).
    % origin_in_nm is 1 x 3, in nm, the corner of the stack in real units
    % spacing_in_nm is 1 x 3, in nm, the spacing between voxels in xyz
    % (the whole level-of-the-octree business is dealt with in client code in this
    % version)
    
    origin_in_um = origin_in_nm / 1e3 ;  %  nm -> um
    spacing_in_um = spacing_in_nm / 1e3 ;  % convert for this level, convert nm -> um
    
    r_in_voxels_shifted = r_in_voxels - 0.5 ;  
        % Not sure where this comes from...
        % My guess is that r_in_voxels contains voxel indices where the
        % index of the lowest-index voxel is (1,1,1) (i.e. Matlab-style), 
        % and the low-coord corner of this voxel is at (0,0,0) in the
        % destination coordinates, so if the spacing vector was [1, 1, 1], then
        % the center of this voxel should be at (0.5, 0.5, 0.5), and this
        % shift makes that so.  That must be it.  --ALT, 2019-06-06
    r_in_um_relative_to_origin = bsxfun(@times, r_in_voxels_shifted, spacing_in_um) ;
    result = bsxfun(@plus, r_in_um_relative_to_origin, origin_in_um) ;
    %result = r_in_voxels .* (ones(size(r_in_voxels,1),1)*spacing_in_um) + ones(size(r_in_voxels,1),1)*origin_in_um ;
    %result_check = origin_in_um + spacing_in_um .* (r_in_voxels - 0.5) ;
    %   result_check always matches result
    %assert(isequal(result, result_check)) ;
end
