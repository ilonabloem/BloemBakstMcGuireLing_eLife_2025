function rootPath = projectRootPath()

% Find the folder that contains this function
filePath    = fileparts(which('projectRootPath'));

% Break up path to find main exp folder
pathComp    = split(filePath, '/');
indx        = strcmp(pathComp, 'Code');

if isempty(pathComp{1}) > 0
    pathComp    = pathComp(2:find(indx)-1);   
else
    pathComp    = pathComp(1:find(indx)-1);
end

% Concatenate path 
rootPath    = filesep;
for n = 1:numel(pathComp)   
    rootPath = strcat(rootPath, sprintf('%s%s', pathComp{n}, filesep));
end

% Add paths and necessary toolboxes
% Dependencies:
% - knkutils
% - circstat-matlab
% - Violinplot-Matlab

if exist(fullfile(rootPath, 'Toolboxes', 'Violinplot-Matlab'), 'dir') > 0
    addpath(genpath(filePath), genpath(fullfile(rootPath, 'Toolboxes')))
else
    addpath(genpath(filePath), genpath(fullfile('/projectnb2', 'vision', 'SpatialUncertainty', 'Toolboxes')))
end