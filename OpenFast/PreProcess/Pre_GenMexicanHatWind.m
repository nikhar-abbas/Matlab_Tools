% Pre_GenMexicanHatWind
%% Generate a mexican hat wind profile. 
% - This runs Pre_GenWindFile.m 


%% Define Inputs
outdir = '/Users/nabbas/Documents/TurbineModels/WindFiles/MexicanHat';
filename = 'MexicanHat_14mps_2A.wnd';

% General wind charactaristics
v_mean = 14;   % mean wind speed
v_max = 2;     % mexican hat peak over mean wind speed
t_shift = 50;  % amount of time to shift full mexican hat profile 


% Mexican hat profile characteristics
f = 1/20; % peak frequency
tlength = 75; % length of signal for mexican hat
dt = 0.01;

% define mexican hat
time = -tlength/2:dt:tlength/2;
mh = (1 - 2*pi^2*f^2.*time.^2).*exp(-pi^2*f^2.*time.^2);

% time shift to start at zero
time = time + tlength/2;

% modify for full time series
t_append = 0:dt:t_shift-dt;
time = [ t_append, time + t_shift ];
mh = [ zeros(1,length(t_append)), mh];

% modify mexican hat to be desired wind profile
windspeeds = v_mean + mh*v_max;
% Generate step wind file
Pre_GenWindFile([outdir filesep filename],time,windspeeds)

