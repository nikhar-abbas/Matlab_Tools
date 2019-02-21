% GenStepwind
%% Generate wind files with this script. 
% - This just runs Pre_GenStepWind.m, but makes setup a bit easier


%% Define Inputs
filename = 'NoShr_5-8_Inc.5_50s.wnd';
windspeeds = [5:.5:8];
tstep = 50;
outdir = pwd;


% Generate wind
Pre_GenStepWind([outdir filesep filename],tstep,windspeeds)

%%
pwd