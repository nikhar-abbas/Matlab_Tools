% GenStepwind
%% Generate wind files with this script. 
% - This runs Pre_GenWindFile.m to generate a step wind file

%% Define Inputs
outdir = '/Users/nabbas/Documents/TurbineModels/WindFiles/StepWind';
filename = 'NoShr_9-14_Inc1_50s.wnd';
windspeeds = [9:14];
tstep = 50;
% tlength = 500;

t = 0;
j = 1;
for i = 1:length(windspeeds)
    time(j) = t;
    windspeedArray(j) = windspeeds(i);
    j = j+1;
    t = t+tstep-1
    
    time(j) = t;
    windspeedArray(j) = windspeeds(i);
    t = t+1;
    j = j+1;
end

% Generate step wind file
Pre_GenWindFile([outdir filesep filename],time,windspeedArray)

