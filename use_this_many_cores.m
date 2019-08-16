function use_this_many_cores(core_count_desired)
    poolobj = gcp('nocreate') ;  % If no pool, do not create new one.
    if isempty(poolobj) ,                
        parpool(core_count_desired) ;        
    else
        core_count = poolobj.NumWorkers ;
        if core_count ~= core_count_desired ,
            delete(poolobj) ;
            parpool(core_count_desired) ;
        end
    end
    poolobj = gcp('nocreate');  % If no pool, do not create new one.  But there should be one at this point...
    core_count = poolobj.NumWorkers ;
    fprintf('Using %d cores.\n', core_count) ;
end
