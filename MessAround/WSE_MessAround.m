%% WSE_MessAround.m
% Lets try to get a wind speed estimator working
%
% The design of this estimator is based off of Torbin Knutsen's Wind Energy
% publication titled "Preditciont models for wind speed at turbine
% locations. Of course, give the man credit, where credit is due.
%
% Nikhar Abbas


%% Load a turbine simulation for some data to work with
% load('/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU10MWRWT/Baseline/Outfiles_Simulink/11.4mps_turb.mat');  % loads simout structure
% load('/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU10MWRWT/Baseline/Outfiles_Simulink/10mps_steady.mat');  % loads simout structure
% simout = simtemp;
% A Cp surface for this turbine
% cpscan = load('/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU10MWRWT/CpScan/MatFiles/CpScan_FullSurf.mat');
cpscan = load('/Users/nabbas/Documents/TurbineModels/NREL_5MW/CpScan/CpScan.mat');

% And some linearization parameters
% addpath(genpath('/Users/nabbas/Documents/TurbineModels/TurbineControllers/SimulinkControllers/TSR_Tracking'))
% addpath('/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU10MWRWT/DTU_10MW_Simulink')
addpath(genpath('/Users/nabbas/Documents/TurbineModels/NREL_5MW/5MW_Land'));
% ContParam = Pre_ContParam_TSR_DTU10MW;
ContParam = Pre_ContParam_TSR_NREL5MW;
[A_v,Bb,GS,Beta_op,vv] = Pre_TSRtracking_GS(ContParam,cpscan);
ContParam.GS = GS;
%% Define turbine parameters and that fun stuff
J = ContParam.J;                    % Rotor Inertia, kg*m^2
rho = ContParam.rho;                % Air Density, kg/m^3
Rr = ContParam.RotorRad;
A = pi*ContParam.RotorRad^2;        % Rotor Swept Area, m^2

% Setup save variables
clear xh xhmat vh nv Vm1
xhmat = [];
vh = [];
Vm1 = zeros(2,2);
%% Load simulation output and Cp data
T = simout.Time;                    % Time, s
% T = T(1:end/4);
Tg = simout.GenTq*1000;             % Generator Torque, N*m
% Load Cp data
TSRvec = cpscan.TSR;
Betavec = cpscan.BlPitch .* pi/180;
Cpmat = cpscan.Cpmat;

% Measured outputs at current time
beta = simout.BldPitch1 .* pi/180; 
tau_g = simout.GenTq .* 1000;
% v_n = simout.RtVAvgxh; % Assume known wind speed, not measured wind speed @ rotor. 
v_n = simout.Wind1VelX;  
om_m = simout.GenSpeed.*(pi/30)./ContParam.GBRatio;
dt = simout.Time(2) - simout.Time(1);

%% The filter
% Initial Conditions
L = 100;
P = diag([0.1, 0.1, 1])^2;
K = zeros(3,1);
H = [1 0 0];
F = [om_m(1), 0, v_n(1)];
v_m = 11.4; v_n(1);
v_t = 0.01;
xh(:,1) = [om_m(1) v_t v_m]';
om_r = om_m(1);

% Covariances
% R = diag([0.02, 1]);
R = 0.02;
V22 = 2^2/600;

% Run the filter
for ti = 1:length(T)

% Find current estimated operating conditions, and saturate
v_r = v_t + v_m;
% v_r = max( min(v_r, vv(end)), vv(1));
TSRe = om_r*Rr/v_r; 
TSRe = max( min(TSRe, TSRvec(end)), TSRvec(1));

% Find current Cp 
CpTSR = zeros(1,length(TSRvec));
for TSRi = 1:length(TSRvec)
    CpTSR(TSRi) = interp1(Betavec, Cpmat(TSRi,:), beta(ti)); % Vector of Cp values corresponding to operational beta
end
Cp = interp1(TSRvec, CpTSR, TSRe)
Cp = max(0,Cp);

% Cp = simout.RtAeroCp(ti);
% Calculate Jacobians
F = WSE_Jacobian(ContParam, Cp, om_r, v_m, v_t, L);



% THE FILTER IN ACTION
% Update Covariance 
tui = 0.1;
V11 = pi * v_m^3 * tui^2 / L;
V = [V11 0; 0 V22];
Q = diag([5e-5, V11, V22]);


%Prediction Update
dxh = state_f(ContParam, tau_g(ti), Cp, v_m, v_r, v_t, om_r, L);
xh = xh + dt*dxh;
% P = F*P*F' + Q; % EKF
dP = F*P + P*F' + Q - K*R*K';
P = P + dt*dP;

% Measurement Update
% y_til = [om_m(ti); v_n(ti)] - [om_r; v_r]
y_til = [om_m(ti)] - [om_r]
y_tilmat(ti) = y_til;
if isnan(y_til)
    break
end
xh = xh + K*y_til;

S = H*P*H' + R;
K = P*H'/R;
% P = (eye(3) - K*H)*P; %EKF
% P = P - K*S*K';
P = (eye(3) - K*H)*P*(eye(3) - K*H)' + K*R*K'; 


% if T < 20
%     v_m = 8;
% end
% Define variables
om_r = xh(1);
v_t = xh(2);
v_m = xh(3);
v_r = v_m + v_t;

% Save 
xhmat = [xhmat xh];
vh = [vh v_r];


end


%% 
function F = WSE_Jacobian(ContParam, Cp, om_r, v_m, v_t, L)
% This finds the jacobian matrix used in the EKF for Wind Speed Estimation
% -This could (should) be made a function that is F(v,omega,beta) without
% the need to keep loading Cp Surface 
%
% Inputs: ContParam - Structure of control parameters
%         Avec - vector system poles, w.r.t. v
%         vv - vector of wind speed linearization points , m/s
%         cpscan - Output from Run_Cpscan.m
%         om_r - rotor speed, rad/s
%         beta - Blade Pitch Angle
%         v_m - estimated mean wind speed - 10 min average
%         v_t - estimated wind speed - turbulent component
%
% Nikhar Abbas

% Load ContParam variables
J = ContParam.J;                    % Rotor Inertia
rho = ContParam.rho;                % Air density
R = ContParam.RotorRad;             % Rotor Radius, R
A = pi*R^2;                         % Rotor Swept Area, m^2
pA = ContParam.GS.pA;

% Slow filter "length"
% L = 300;

% F = [-0.08, 1/(2*J)*rho*A*Cp*3*v_t^2*1/om_r, 1/(2*J)*rho*A*Cp*3*v_m^2*1/om_r;...
F = [pA(1)*(v_m+v_t) + pA(2), 1/(2*J)*rho*A*Cp*3*(v_m+v_t)^2*1/om_r, 1/(2*J)*rho*A*Cp*3*(v_m+v_t)^2*1/om_r;...   
    0, pi*v_m/(2*L), pi*v_t/(2*L);...
    0, 0, 0];
    



end

%%
function [xhd, tau_r] = state_f(ContParam, tau_g, Cp, v_m, v_r, v_t, om_r, L)
% Updated the state estimate
% - This needs more commenting
% Nikhar Abbas

% Load ContParam variables
J = ContParam.J;                    % Rotor Inertia
rho = ContParam.rho;                % Air density
R = ContParam.RotorRad;             % Rotor Radius, R
A = pi*R^2;                         % Rotor Swept Area, m^2
Ng = ContParam.GBRatio;             % Gearbox Ratio
% Slow filter "length"
% L = 300;

tsr = om_r*ContParam.RotorRad/v_r;

% Torque calculation
tau_r = 1/2*rho*A*R*Cp*v_r^2 * 1/tsr;
% Wind variance
a = pi*v_m/(2*L);

% generate noise
% nv = V* randn(2,1);

% state
omd = 1/J * (tau_r - Ng*tau_g);
vtd = -a*v_t;
vmd = 0;

xhd = [omd vtd vmd]';






end