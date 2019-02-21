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
Cpl = linavg.C(25,:); % Generator Speed as output (rad/s)
Dpl = linavg.D(25,1);

G = tf(ss(Apl,Bpl,Cpl,Dpl));

% G = tf([-2.982e-6],[1 0.106]);


%% Make a controller 
% -- Parameters for DTU 10MW --
J = 156348032.208;                  % Rotor Inertia (kg-m^2)
rho = 1.225;                        % Air Density (kg/m^3)
% operating conditions
v = 6;
tsr = 7.5; 10.85;
Cp = 0.483;
R = 89.15;                          % Rotor Radius (m)
Ar = pi*R^2; 

% Parameters
A = 1/(2*J) * rho*Ar*R^2*v*Cp/(tsr^2);

% d(\dpt{omega})
domd = 1/J * .5*rho*Ar*R^2*v*(1/tsr^2)* (.0555 - 2*Cp) % .0555 is estimation of dCp/dlambda


% A = -0.106;
B = -50/J * 50 * 30/pi ;

% b = Bd*Cd;


% ---- Controller ----
% Calc gains
zeta = .7;                   % more than critically damped - try to account for zero in CL TF
% Ts = 2;                    % settling time (s)
Tr = 15;
% om_n = 4.6/(zeta*Ts);       % Natural Frequency
om_n = 1.8/Tr;

% VS_Kp = (2*zeta*om_n - Ad)/(Bd*Cd)
% VS_Ki = om_n^2/(Bd*Cd)

Kp = 1/B * (2*zeta*om_n + A);
Ki = om_n^2/B;

s = tf('s');
H = Kp + Ki/s;


sys = minreal(feedback(H*G,1))

%% Some plots

step(sys)