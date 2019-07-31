% ControlEnergy
% This script is written to do some analysis on the control energy
% necessary to move a system from one state to another. We use the
% controllability gramian and a define domain to do this analysis. For
% finite time analysis, we include an argument for time

%% Save figures?
savefigs = 1
spath = '/Users/nabbas/Documents/Publications/WindTech2019/Controllability_Manuscript/Figures/';

%% Find Control Energy
thetavec = 0;
% thetavec = [35, 45, 55];
for thi = 1:length(thetavec)
%% Define System - MassDamper
theta = thetavec(thi); % can't be 90, Gramian can't be calculated

m = 10;
k1 = 5;
k2 = k1;
c1 = 2;
c2 = c1;
A = [0 1 0 0;...
    -k1/m -c1/m 0 0;...
    0 0 0 1;...
    0 0 -k2/m -c2/m];
B = [0 0;...
    1/m sind(theta)/m;...
    0 0;...
    0 cosd(theta)/m];
C = eye(4);
D = zeros(4,2);

sys = ss(A,B,C,D);

%% Start t-loop
tvec = [18 20 22];
% tvec = 20;
for ti = 1:length(tvec)
    
%% Parameters for this run
t = tvec(ti);
%% Controllability Gramian, and SVD
opt = gramOptions('TimeIntervals',[0 t]);
Wc = gram(sys,'c',opt);
[U,S,V] = svd(Wc);

%% Find Energy Surface - MassDamper
xvec = linspace(-1,1,31); % need to be odd
yvec = linspace(-1,1,31);

x0 = [1 0 2 0]';
xdes = [0 0 0 0]';

Esurf = EnergySurface_Data(Wc,A,t,x0,xdes,xvec,yvec);

%% Make Plots

figure(11)
sp = subplot(1,3,ti);
% sp = subplot(1,3,thi);
SurfPlot(xvec, yvec, Esurf');
caxis([0 0.5])
% caxis([0 2])
sppos = get(sp,'Position');
if ti == length(tvec) && thi == length(thetavec)
    colorbar('EastOutside')
end
tstring = [num2str(t), ' seconds'];
% tstring = [num2str(theta), ' Degrees'];
title(tstring, 'fontsize', 9, 'FontWeight', 'normal')
sp.Position
sp.Position(3) = sppos(3);
sp.OuterPosition(4) = sppos(4) + 0.07;
sp.Position
if savefigs
%     saveas(gcf, [spath 'MassDamper_tfix.eps'], 'epsc2')
    saveas(gcf, [spath 'MassDamper_thetafix.eps'], 'epsc2')
end

sgtitle('Control Energy to Return to the Origin (Nm)', 'FontName', 'Times');

%% End t-loop
end

%% End theta-loop
end









%% ----------------
% Floating Platform analysis
%% ----------------

%% Find Control Energy
Lindir = '/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU10MWRWT_NAUTILUS_GoM_FAST_v1.00/Linearizations/Case1_AllFloat';
LinfileBase = 'DTU_10MW_NAUTILUS_GoM'
nlin = 24;

% Load cases
% [linout, linavg] = Post_GetLinear(Lindir, LinfileBase, nlin)
[linout, linavg] = Post_LinAnalysis_1(Lindir, LinfileBase, nlin);
% Load system
A = linavg.A;
B = linavg.B;
C = linavg.C;
D = linavg.D;

sys = ss(A,B,C,D);

%% Start t-loop
% tvec = [1.58 3.16 4.74];
p1 = 1.58;
tvec = [p1*2 p1*3 p1*4];
% tvec = 20;
for ti = 1:length(tvec)
    
%% Parameters for this run
t = tvec(ti);
%% Controllability Gramian, and SVD
opt = gramOptions('TimeIntervals',[0 t]);
Wc = gram(sys,'c',opt);
[U,S,V] = svd(Wc);

%% Find Energy Surface - MassDamper
xvec = linspace(-2,2,31); % need to be odd
yvec = linspace(-3,3,31);

x0 = zeros(4,1);
x0(1) = 1;
x0(2) = 2;
xdes = zeros(4,1);

Esurf = EnergySurface_Data(Wc,A,t,x0,xdes,xvec,yvec);

%% Make Plots

figure(21)
sp = subplot(1,3,ti);
SurfPlot(xvec, yvec, Esurf');
xlabel('Platform Roll Angle (deg)')
ylabel('Platform Pitch Angle (deg)')

caxis([1e5 .5e9])
% caxis([0 7e4])
sppos = get(sp,'Position');
if ti == length(tvec)
    colorbar('EastOutside')
end
tstring = [num2str(t), ' seconds'];
title(tstring, 'fontsize', 9, 'FontWeight', 'normal')
sp.Position
sp.Position(3) = sppos(3);
sp.OuterPosition(4) = sppos(4) + 0.07;
sp.Position
if savefigs
    saveas(gcf, [spath 'Nautilus_ContE.eps'], 'epsc2')
end

sgtitle('Control Energy to Return to the Vertical', 'FontName', 'Times');

%% End t-loop
end




%% Compare Nautilius and OO Semi-submersibles
A1 = A;
sys1 = sys;         % Nautilus - need to run previous section

% Load OO-semi
Lindir = '/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU10MWRWT_OO_GoM_FAST_v1.00/Linearizations/Case1_AllFloat';
LinfileBase = 'DTU_10MW_OO_GoM'
nlin = 24;

% Load cases
% [linout, linavg] = Post_GetLinear(Lindir, LinfileBase, nlin)
[linout, linavg] = Post_LinAnalysis_1(Lindir, LinfileBase, nlin);

A = linavg.A; A2 = A;
B = linavg.B;
C = linavg.C;
D = linavg.D;
sys2 = ss(A,B,C,D);

%% Controllability Gramians, and Energy
t = p1*3;
opt = gramOptions('TimeIntervals',[0 t]);
Wc1 = gram(sys1,'c',opt);
Wc2 = gram(sys2,'c',opt);

Esurfm(:,:,1) = EnergySurface_Data(Wc1,A1,t,x0,xdes,xvec,yvec);
Esurfm(:,:,2) = EnergySurface_Data(Wc2,A2,t,x0,xdes,xvec,yvec);

%% Make Plot
figure(31)
enum = 2;
for esi = 1:enum
sp = subplot(1,enum,esi);
SurfPlot(xvec, yvec, Esurfm(:,:,esi)');
xlabel('Platform Roll Angle (deg)')
ylabel('Platform Pitch Angle (deg)')

caxis([1e5 2e8])
sppos = get(sp,'Position');
if esi == enum
    colorbar('EastOutside')
end
tstring = {'NAUTILUS-10','OO-Star'};
title(tstring{esi}, 'fontsize', 9, 'FontWeight', 'normal')
sp.Position
sp.Position(3) = sppos(3);
sp.OuterPosition(4) = sppos(4) + 0.07;
sp.Position

sgtitle('Control Energy to Return to Vertical', 'FontName', 'Times');

set(gcf,'PaperPosition',[0 0 5 2])

if savefigs  
    saveas(gcf, [spath 'PtfmComp_ContE.eps'], 'epsc2')
end
end





%% SVD Plot
Lindir = '/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU10MWRWT_NAUTILUS_GoM_FAST_v1.00/Linearizations/Case8_PitSur';
LinfileBase = 'DTU_10MW_NAUTILUS_GoM'
nlin = 24;
[linout, linavg] = Post_LinAnalysis_1(Lindir, LinfileBase, nlin);
% Load system
A = linavg.A;
B = linavg.B;
C = linavg.C;
D = linavg.D;

sys_svd = ss(A,B,C,D);

% Gramian and SVD
opt = gramOptions('TimeIntervals',[0 p1*10]);
Wc_SVD = gram(sys_svd,'c');
[U,S,V] = svd(Wc_SVD);


figure(41)
plot(diag(S),'x-','linewidth',1.5)
xlim([0 length(A)])
xlabel('Singular Value Index')
ylabel('Singular Value')
title('DTU 10MW of NAUTILUS-10 Singular Values')
grid on
set(gca,'FontSize',8,'FontName','Times')
set(gcf,'PaperPosition',[0 0 4.5 2])

if savefigs
    saveas(gcf, [spath 'SingularValues.eps'], 'epsc2')
end


%%
function Esurf = EnergySurface_Data(Wc,A,t,x0,xdes,xvec,yvec)
% Find a surface of energy to return the origin within some amount of time.
% 
% Inputs:
%   Wc - Controllability Gramian
%   A - System matrix
%   t - time
%   x0 - initial state. Input value of 1 should denote x value, 2 denotes y
%        value
%   xdes - final state
%   xvec - vector of xvalues denoting surface
%   yvec - vector of yvalues denoting surface

%% Make a surface of Energy to get "Home"
    xs = length(xvec);
    ys = length(yvec);
    
    % Find x,y, variables
    xloc = find(x0 == 1);
    yloc = find(x0 == 2);
    
    
    for xi = 1:xs
        for yi = 1:ys
            x0(xloc) = xvec(xi);
            x0(yloc) = yvec(yi);
%              x0 = [xloc 0 yloc 0]';           % Starting location
%              xdes = [0 0 0 0]';                 % Final Location

            E = (expm(A.*t)*x0 - xdes)'/Wc*(expm(A.*t)*x0 - xdes);

            Esurf(xi,yi) = E;
        end
    end
end

%% make surface plots
function SurfPlot(x,y,z)

s = pcolor(x,y,z);
% s.FaceColor = 'interp';
view(2)
% colorbar
% colormap('hsv')
% xlabel('x-displacement (m)')
% ylabel('y-displacement (m)')

set(gca,'FontSize',8,'FontName','Times')
set(gcf,'PaperPosition',[0 0 7.5 2.2])

end
