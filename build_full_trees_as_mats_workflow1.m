function build_full_trees_as_mats_workflow1(configuration_file_path)

    options = configparser(configuration_file_path);
    % if ~isfield(options,'sampling')
    %     options.sampling = 'uni';
    % end
    myh5 = options.inputh5 ;
    myh5prob = options.h5prob ;

    [brainSize,~,~,~] = h5parser(myh5,myh5prob);

    origin = h5read(options.inputh5,[options.h5prob,'_props/origin']);
    spacing = h5read(options.inputh5,[options.h5prob,'_props/spacing']);
    level = h5read(options.inputh5,[options.h5prob,'_props/level']);
    params.outsiz = brainSize;
    params.ox = origin(1);
    params.oy = origin(2);
    params.oz = origin(3);
    params.sx = spacing(1);
    params.sy = spacing(2);
    params.sz = spacing(3);
    params.level = level;

    params.voxres = [params.sx params.sy params.sz]/2^(params.level)/1e3; % in um
    options.params = params;

    [subs,~,A,~] = skel2graph(options) ;

    Gin = graph(max(A,A')) ;
    workflow1_full_trees_only_as_mats(Gin,subs,options) ;
end
