function [swcout,swcfiles] = loadSwcFolder(swcfolder)
%%
% loads all swcs in a folder
swcfiles = dir(fullfile(swcfolder,'*.swc'));
swcout = cell(1,length(swcfiles));
parfor iswc = 1:length(swcfiles)
    [swcData,offset,~, ~] = loadSWC(fullfile(swcfolder,swcfiles(iswc).name));
    swcData(:,3:5) = swcData(:,3:5) + ones(size(swcData,1),1)*offset;
    swcout{iswc} = swcData;
end
