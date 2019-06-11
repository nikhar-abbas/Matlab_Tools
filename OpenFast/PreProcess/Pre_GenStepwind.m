% GenStepwind
%% Generate wind files with this script. 
% - This just runs Pre_GenStepWind.m, but makes setup a bit easier


%% Define Inputs
outdir = '/Users/nabbas/Documents/TurbineModels/WindFiles/StepWind';
filename = 'NoShr_9-14_Inc1_50s.wnd';
windspeeds = [9:14];
tstep = 50;


% Generate wind
Pre_GenStepWindFile([outdir filesep filename],tstep,windspeeds)

