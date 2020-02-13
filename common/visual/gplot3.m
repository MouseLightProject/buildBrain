function h = gplot3(varargin)
%function h = gplot3(A, xyz, varargin)
    %GPLOT Plot graph (nodes and edges).
    %   GPLOT(A, xyz) plots the graph specified by the adjacency matrix,
    %   A, and the n-by-3 coordinate array, xyz.
    %   
    %   GPLOT(A, xyz, linespec) uses line type and color specified in the
    %   string LineSpec. See PLOT for possibilities.
    %
    %   h = GPLOT(A, xyz) returns the a handle to the graph.
    %   
    %   h = GPLOT(A, xyz, 'LineWidth', 5, ...) also takes arbitrary arguments
    %   for line properties
    
    remaining_args = varargin ;
    
    if length(remaining_args)>=1 ,
        arg = remaining_args{1} ;
        if isgraphics(arg, 'axes') ,
            ax = arg ;
            remaining_args = remaining_args(2:end) ;
        else
            ax = gca() ;
        end
    end       
    
    if length(remaining_args)>=1 ,
        A = remaining_args{1} ;
        remaining_args = remaining_args(2:end) ;
    else
        error('Need an adjacency matrix') ;
    end

    if length(remaining_args)>=1 ,
        xyz = remaining_args{1} ;
        remaining_args = remaining_args(2:end) ;
    else
        error('Need coordinate matrix') ;
    end
    
    % Returns i and j, lists of connected nodes
    [i,j] = find(A);

    % Extact 
    X = [ xyz(i,1) xyz(j,1)]';
    Y = [ xyz(i,2) xyz(j,2)]';
    Z = [ xyz(i,3) xyz(j,3)]';
    
    % Add NaN values to break between line segments
    X = [X; NaN(size(i))'];
    Y = [Y; NaN(size(i))'];
    Z = [Z; NaN(size(i))'];

    % Serialize the x and y data
    X = X(:);
    Y = Y(:);
    Z = Z(:);
    
    % If linespec given, then use it
    if mod(length(remaining_args),2) == 1 ,
        % There's an odd number of remaining args, so assume the first is a
        % linespec
        linespec = remaining_args{1} ;
        remaining_args = remaining_args(2:end) ;
        h = plot3(ax, X, Y, Z, linespec);        
    else
        % There's a an even number of remaining args, so plot without
        % linespec
        h = plot3(ax, X, Y, Z);
    end
        
    % Now apply the rest of the args
    for i = 1:2:length(remaining_args) ,
        set(h, remaining_args{i}, remaining_args{i+1});
    end    
end
