function eout = buildgraph(Asub, source)
    % Asub should be the adjacency matrix for an undirected graph, and thus
    % symmetric.
    if nargin<2 ,
        source = find(sum(Asub)==1, 1) ;
    end
    [~,~,pred] = graphshortestpath(Asub, source, 'directed', false) ;
    eout = [((1:size(Asub,1))') (pred(:))];
    eout(eout(:,2)==0,:) = [];
end