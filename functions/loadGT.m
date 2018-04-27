function [swcpixlocs,conn,swcfilepath] = loadGT(swcfile,params,index)
%% GT DATA
SX = params.sx;
SY = params.sy;
SZ = params.sz;
voxres = [SX SY SZ]/2^(params.level)/1e3; % in um
params.voxres = voxres;
numNeuron = length(swcfile);
if nargin<3
    index = 1:numNeuron;
end
if 1
    clear swcpixlocs
    numNeuron = length(swcfile);
    swcpixlocs = cell(1,numNeuron);
    swcfilepath = cell(1,numNeuron);
    parfor ii=1:numNeuron
        if isempty(swcfile{ii})
            continue
        end
        %%
        sprintf('Loading: %s',swcfile{ii})
        [swcData,offset,color, header] = loadSWC(swcfile{ii});
        swcData(:,3:5) = swcData(:,3:5) + ones(size(swcData,1),1)*offset;
        swclocs = swcData(:,3:5)*1000;
        swclocs = swclocs-ones(size(swclocs,1),1)*[params.ox params.oy params.oz];
        swcpixlocs{ii} = round(swclocs./(ones(size(swclocs,1),1)*(voxres*2^0))/1000);
        %%
        % edges
        edges = swcData(:,[1 7]);
        edges(any(edges==-1,2),:) = [];
        if isempty(edges)
            E=[];
        else
            E = sparse(edges(:,1),edges(:,2),1,max(edges(:)),max(edges(:)));
        end
        conn{ii} = max(E,E');
        swcfilepath{ii} = swcfile{ii};
    end
%     numNeuron = length(swcpixlocs);
%     somalocs = zeros(numNeuron,3);
%     for ii=1:numNeuron
%         if ~isempty(swcpixlocs{ii})
%             somalocs((ii),:) = swcpixlocs{ii}(1,:);
%         end
%     end
else
    %%
    % %% somalocs are from worksheet
    % somalocs=[75697.3, 42460.0, 19956.5
    %     76267.1, 42106.1, 20354.4
    %     76323.5, 41874.5, 20146.9
    %     75971.7, 42292.1, 20006.6
    %     76463.8, 42439.8, 19799.1
    %     75602.1, 42240.6, 18966.3
    %     75598.2, 42052.6, 19757.0
    %     76862.7, 41725.7, 19855.2
    %     75530.8, 41984.1, 20217.1
    %     76243.6, 41781.7, 20243.1
    %     76548.9, 42253.8, 19193.8
    %     76492.8, 42688.4, 19609.7
    %     78934.0, 42506.1, 18251.7
    %     79175.1, 42838.9, 18480.2
    %     79140.3, 42853.2, 18187.5
    %     79393.6, 42706.3, 18363.9
    %     78932.8, 42504.9, 18243.7
    %     79025.9, 42206.6, 18190.0
    %     78749.2, 42184.9, 18055.2
    %     79275.3, 42594.7, 17991.1];
    % somalocs = somalocs-ones(size(somalocs,1),1)*[params.ox params.oy params.oz]/1000;
    % somalocs = round(somalocs./(ones(size(somalocs,1),1)*(voxres*2^0)));
end