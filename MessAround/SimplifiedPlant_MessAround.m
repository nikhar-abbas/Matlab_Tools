Run%% Simplified plant model mess-around

%% Lets load some linearization data
% Linearizatin data 
% Lindir = '/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU10MWRWT/Linearizations/AboveRated/Case1';
% LinfileBase = 'DTU_10MW_RWT';

Lindir = '/Users/nabbas/Documents/TurbineModels/NREL_5MW/5MW_Land/Linearization/';
LinfileBase = '5MW_Land';
nlin = 24;
[linout, linavg] = Post_GetLinear(Lindir, LinfileBase, nlin);


%% Plant Model
Apl = linavg.A;
Bpl_t = linavg.B(:,8); % Generator Torque as input
Bpl_b = linavg.B(:,9)*pi/180; % Blade Pitch as input

% Cpl = linavg.C(20,:); % Rotor Speed as output (rpm)
Cpl = linavg.C(25,:)*pi/30; % Generator Speed as output (rad/s)
Ctsr = linavg.C(44,:); 
Dpl = linavg.D(44,8);

G = tf(ss(Apl,Bpl_t,Cpl,Dpl));
G_b = tf(ss(Apl,Bpl_b,Cpl,Dpl));

% 11.5 - step = -8.06e-4
% 13 - step = -0.544
% 20 - step = -0.0129
%% Make a controller 
% -- Parameters for DTU 10MW --
Jr = 156348032.108;                  % Rotor Inertia (kg-m^2)
Jg = 1500.5;
J = Jr + Jg*50^2;
rho = 1.225;                        % Air Density (kg/m^3)
R = 89.05;                          % Rotor Radius (m)
Ar = pi*R^2; 
Ng = 50; 

% operating conditions
v = 11.4;
Bv = 0:20;
% for lind = 1:length(Bv)
tsr = 1.005*R/v; 
Beta = 0;

%%% THIS ISN'T REALLY WORKING RIGHT NOW 4/30 - but the real stuff does ;)
% Cp data
cpscan = load('/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU10MWRWT/CpScan/MatFiles/CpScan_FullSurf.mat');
TSRvec = cpscan.TSR;
Cpvec = cpscan.Cpmat(:,(cpscan.BlPitch == 0))';
Betavec = cpscan.BlPitch;
Cpmat = cpscan.Cpmat;

% Define Cp Operating Points
Vci = 4;
Vco = 24;
RRspeed = 9.6*pi/30;

TSRr = RRspeed*R/11.4;
vv_br = [4:.2:Vci-eps]; 
vv_ar = [Vci:.2:Vco];
vv_ar = vv_ar(2:end);
vv = [vv_br vv_ar];
TSR_br = ones(1,length(vv_br)) * TSRvec(Cpvec == max(Cpvec));
TSR_ar = RRspeed.*R./vv_ar;
TSRop = [TSR_br TSR_ar];
Cpr = interp1(TSRvec,Cpvec,TSR_ar(1));
Cp_op_br = ones(1,length(vv_br)) * max(Cpvec);
Cp_op_ar = Cpr.*(TSR_ar./TSRr).^3;

Cp_op = [Cp_op_br Cp_op_ar];

for toi = 1:length(TSRop)
tsr = TSRop(toi);

% Saturate initial TSR and Beta
Beta = max(min(Beta,Betavec(end)),Betavec(1));
tsr = max(min(tsr,TSRvec(end)),TSRvec(1));

% ---- Cp Operating conditions ----
CpTSR = zeros(1,length(Betavec));
CpB = zeros(1,length(TSRvec));
for Bi = 1:length(Betavec)
    CpTSR(Bi) = interp1(TSRvec, Cpmat(:,Bi), tsr); % Vector of Cp values corresponding to operational TSR
end
for TSRi = 1:length(TSRvec)
    CpB(TSRi) = interp1(Betavec,Cpmat(TSRi,:),Beta); % Vector of Cp values corresponding to operational Beta
end

% Difference vectors
dB = Betavec(1:end-1) + (Betavec(2) - Betavec(1))/2;
dTSR = TSRvec(1:end-1) + (TSRvec(2) - TSRvec(1))/2;
dCp_B = diff(CpTSR)./diff(Betavec); % Difference of Cp w.r.t. Beta
dCp_tsr = diff(CpB)./diff(TSRvec); % Difference of Cp w.r.t. TSR

% Saturate TSR and Beta
Beta = max(min(Beta,dB(end)),dB(1));
tsr = max(min(tsr,dTSR(end)),dTSR(1));

% Find operational parameters
Cp = interp1(Betavec,CpTSR,Beta);
dCpdB = interp1(dB,dCp_B,Beta);
dCpdTSR = interp1(dTSR,dCp_tsr,tsr);

% CpM = max(CpTSR);
% CpTSR2 = CpTSR;
% for CtTi = 1:length(CpTSR)
%     if CpTSR(CtTi) > Cp_op(toi)
% %         CpTSR2(CtTi) = CpM + (CpM - Cp_op(toi))/CtTi
%             CpTSR2(CtTi) = CpM + (CpM - Cp_op)/((Betavec(2)-Betavec(1))*CtTi);
%     end
%     if CpTSR(CtTi) == CpM
%         break
%     end
% end

% That new-new
Beta_op(toi) = interp1(CpTSR,Betavec,Cp_op(toi));
dCpdBvec(toi) = interp1(dB,dCp_B,Beta_op(toi));

dtdb(toi) = 0.5*rho*Ar*R*dCpdBvec(toi)*vv(toi)^3 * 1/(9.6*pi/30);


%
dtdl = 1/(2)*rho*Ar*R*v^2*(1/tsr^2)* (dCpdTSR*tsr - Cp); % assume operating at optimal
dldo = R/v;
dtdo = dtdl*dldo;

A = dtdo/J;
B_t = -50/J ;
B_b = 1/(2*J)*rho*Ar*R*v^2*(1/tsr)*dCpdB * Ng; 

% Wind Disturbance Input gain
dldv = -tsr/v; 
dtdv = dtdl*dldv;
B_v = dtdv/J;

% . ---  controller tuning ---
% Calc gains
zeta = 0.7; % settling time (s)
% Ts = 30; %2;
% Tr = 30;
% om_n = 4.6/(zeta*Ts);       % Natural Frequency
% om_n = 1.8/Tr;
om_n = 0.7;

Kp_b = 1/B_b * (2*zeta*om_n + A);
Ki_b = om_n^2/B_b;

% Kp(Bind) = Kp_b;
% Ki(Bind) = Ki_b;
% end

% Kp = 1/(Bpl*Cpl) * (2*zeta*om_n + Apl);
% Ki = om_n^2/(Bpl*Cpl);
s = tf('s');
H_b = Kp_b + Ki_b/s;

% myplot(Betavec,CpTSR);
% hold on
% myplot(interp1(CpTSR,Betavec,Cp_op(toi)),Cp_op(toi),'x','linewidth',8)
% sys = zpk(minreal((H_b*G_b)/(1 + H_b*G_b)));
end
%% Some plots

% step(sys)
cps = cpscan.Cpmat;
cps(cps<0)=0;
surface(cpscan.BlPitch,cpscan.TSR,cps), hold on
plot3(Beta_op*180/pi, TSRop,Cp_op+0.1,'r','linewidth',2.5)


%% Finding poles and gain schedule

[Avec,Bbvec,Beta_op,vv] = Pre_TSRtracking_Plant(ContParam,cpscan);

% Trim for BldPitch Controller
% Bopind = find(Beta_op>0);
% Avec_BPC = Avec(Bopind(1)-1:end);
% Bbvec_BPC = Bbvec(Bopind(1)-1:end);
% Betaop_BPC = Beta_op(Bopind(1)-1:end);
% vv_bpc = vv(Bopind(1)-1:end);
Bopind = find(Beta_op>0);
Avec_BPC = Avec(Bopind(1):end);
Bbvec_BPC = Bbvec(Bopind(1):end);
Betaop_BPC = Beta_op(Bopind(1):end);
vv_bpc = vv(Bopind(1):end);



% Blade Pitch gains
% Desired behavior
PC_zeta = ContParam.PC_zeta;
PC_om_n = ContParam.PC_om_n;

% Calculate Gains
Kp = 1./Bbvec_BPC .* (2*PC_zeta*PC_om_n + Avec_BPC);
Ki = PC_om_n^2./Bbvec_BPC ;