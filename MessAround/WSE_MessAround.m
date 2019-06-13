%% WSE_MessAround.m
% Lets try to get a wind speed estimator working
%
% The design of this estimator is based off of Torbin Knutsen's Wind Energy
% publication titled "Preditciont models for wind speed at turbine
% locations. Of course, give the man credit, where credit is due.
%
% Nikhar Abbas


%% Load a turbine simulation for some data to work with
load('/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU10MWRWT/Baseline/Outfiles_Simulink/11.4mps_turb.mat');  % loads simout structure

% And a related cp surface
cpscan = load('/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU10MWRWT/CpScan/MatFiles/CpScan_FullSurf.mat');

%% Define turbine parameters and that fun stuff
J = ContParam.J;                    % Rotor Inertia, kg*m^2
rho = ContParam.rho;                % Air Density, kg/m^3
A = pi*ContParam.RotorRad^2;        % Rotor Swept Area, m^2

%% Load simulation output data
t = simout.Time;                    % Time, s
Tg = simout.GenTq*1000;             % Generator Torque, N*m




















function Tr = RotTorque(lambda,beta,om_m)






end