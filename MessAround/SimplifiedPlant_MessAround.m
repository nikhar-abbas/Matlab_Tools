%% Simplified plant model mess-around

%% Lets load some linearization data
% Linearizatin data 
Lindir = '/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU10MWRWT/Linearizations/BelowRated/Case1';
LinfileBase = 'DTU_10MW_RWT';
nlin = 24;
[linout, linavg] = Post_GetLinear(Lindir, LinfileBase, nlin);


%% Plant Model
Apl = linavg.A;
Bpl = linavg.B(:,8); % Generator Torque as input
% Cpl = linavg.C(20,:); % Rotor Speed as output (rpm)
Cpl = linavg.C(25,:)*pi/30; % Generator Speed as output (rad/s)
Dpl = linavg.D(25,1);

G = tf(ss(Apl,Bpl,Cpl,Dpl));


%% Make a controller 
% -- Parameters for DTU 10MW --
Jr = 156348032.108;                  % Rotor Inertia (kg-m^2)
Jg = 1500.5;
J = Jr + Jg*50^2;
rho = 1.225;                        % Air Density (kg/m^3)
% operating conditions
v = 6;
tsr = 7.5; 7.5; 10.85;
Cp = 0.47; 8;
R = 89.05;                          % Rotor Radius (m)
Ar = pi*R^2; 

% Cp data
cpscan = load('/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU10MWRWT/CpScan/MatFiles/CpScan_0bl_HiRes.mat');
TSRvec = cpscan.TSR;
Cpvec = cpscan.Cpmat';

% Find Cp Operating conditions 
dCp = diff(Cpvec)./diff(TSRvec);
dTSR = TSRvec(1:end-1) + (TSRvec(2) - TSRvec(1))/2;
if tsr > dTSR(end)
    tsr = dTSR(end);
elseif tsr < dTSR(1)
    tsr = dTSR(1);
end

Cp = interp1(TSRvec,Cpvec,tsr);
dCpdTSR = interp1(dTSR,dCp,tsr);

% . ---  controller tuning ---
dtdl = 1/(2)*rho*Ar*R^2*v*(1/tsr^2)* (dCpdTSR*tsr - Cp); % assume operating at optimal
dldo = R/v;
dtdo = dtdl*dldo;

A = dtdo/J;

B = -50/J ;

% --- wind speed response (disturbance) ---
dldv = -tsr/v; 
dtdv = dtdl*dldv;
Bd = dtdv/J;

% ---- Controller ----
% Calc gains
zeta = 0.7; % settling time (s)
% Ts = 30; %2;
% Tr = 30;
% om_n = 4.6/(zeta*Ts);       % Natural Frequency
% om_n = 1.8/Tr;
om_n = 1/20;

Kp = 1/B * (2*zeta*om_n + A);
Ki = om_n^2/B;

% Kp = 1/(Bpl*Cpl) * (2*zeta*om_n + Apl);
% Ki = om_n^2/(Bpl*Cpl);
s = tf('s');
H = Kp + Ki/s;


sys = zpk(minreal((H*G)/(1 + H*G)));

%% Some plots

step(sys)


%% Possible Poles


Bp = -3.123e-07;
Ap = 0.03791;
plant = zpk( (Bp*(Kp*s + Ki))/(s^2 + (Ap + Kp*Bp)*s + Bp*Ki) )


