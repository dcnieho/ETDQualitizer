function createETDQ_MLTBX(prjFile, toolboxVersion)
% Package toolbox as MLTBX file.
% createToolboxMLTBX(prjFile, toolboxVersion) builds the MLTBX file and saves
% it in the release folder. Input prjFile is the name of the toolbox
% packaging file and toolboxVersion is a string of the form Major.Minor.Bug.Build.

if ~isfile(prjFile)
    error("Unable to find " + "'" + prjFile+ "'");
end

% Build file tree
p = fileparts(mfilename("fullpath"));
copyfile(fullfile(p,'..','matlab'),fullfile(p,'matlab'))
mkdir(fullfile(p,'matlab','examples'))
copyfile(fullfile(p,'..','example','matlab.m'),fullfile(p,'matlab','examples','ETDQualitizer_demo.m'))
copyfile(fullfile(p,'demos.xml'),fullfile(p,'matlab','demos.xml'))
copyfile(fullfile(p,'..','example','data'),fullfile(p,'matlab','examples','data'))

addpath(genpath(fullfile(p,'matlab')))
publish(fullfile(p,'matlab','examples','ETDQualitizer_demo.m'))

% Load packaging info
packagingData = matlab.addons.toolbox.ToolboxOptions(prjFile);

% Update the version number
if ~strcmp(toolboxVersion,'master') % if triggered by workflow distpatch, ignore
    packagingData.ToolboxVersion = toolboxVersion;
end

% Set name of output file to that required by FEX
packagingData.OutputFile = strcat(toolboxVersion, ".mltbx");

% Create toolbox MLTBX
matlab.addons.toolbox.packageToolbox(packagingData);
