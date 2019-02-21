%% Spring Mass Damper
% Lets play with a spring-mass-damper to see what we can learn about basic
% physics

%% Setup Parameters
m = 5;              % mass
k = 8;              % spring constant
b = 20;              % damping constant


%% make a state space representation
A = [0 1; -b/m -k/m];
B = [0; 1/m];
C = eye(2);
D = [0;0];

smd = ss(A,B,C,D);

B2 = [0 .2; 1/m 0];
D2 = [0 0; 0 0];

smd2 = ss(A,B2,C,D2);

%% Find controllability Grammian and related SVD
Wc = gram(smd,'c');
[U,S,V] = svd(Wc);
rot = V;


%% Plot ellipsoid - in steps
close all
theta = atan2(rot(2,1),rot(1,1));
% PlotEllipse(S(1,1),S(2,2),theta,Wc);

I = eye(2);
VI = V'*I;
SVI = S*VI;
USVI = U*SVI;
US = U*S;

figure(1)
% subplot(2,2,1)
myplot([0 I(1,1)], [0 I(2,1)]); hold on
myplot([0 I(1,2)], [0 I(2,2)]);
PlotEllipse(1,1)
xlabel('$x$','Interpreter','Latex')
ylabel('$\dot{x}$','Interpreter','Latex')
axis equal

figure
% subplot(2,2,3)
myplot([0 VI(1,1)], [0 VI(2,1)]); hold on
myplot([0 VI(1,2)], [0 VI(2,2)]);
PlotEllipse(1,1)
xlabel('$x$','Interpreter','Latex')
ylabel('$\dot{x}$','Interpreter','Latex')
axis equal

figure
% subplot(2,2,4)
myplot([0 SVI(1,1)], [0 SVI(2,1)]); hold on
myplot([0 SVI(1,2)], [0 SVI(2,2)]);
PlotEllipse(S(1,1),S(2,2))
myplot([0 S(1,1)], [0, S(2,1)],'k--','linewidth',1);
text(S(1,1)/2, S(2,1)/2 - .002, '\sigma_1', 'fontsize',18)
myplot([0 S(1,2)], [0, S(2,2)],'k--','linewidth',1);
text(S(1,2)/2 + .001, S(2,2)/2, '\sigma_2', 'fontsize',18)
xlabel('$x$','Interpreter','Latex')
ylabel('$\dot{x}$','Interpreter','Latex')
axis equal

figure
% subplot(2,2,2)
myplot([0 USVI(1,1)], [0 USVI(2,1)]); hold on
myplot([0 USVI(1,2)], [0 USVI(2,2)]);
PlotEllipse(S(1,1),S(2,2),theta)
myplot([0 US(1,1)], [0, US(2,1)],'k--','linewidth',1);
text(US(1,1)/2, US(2,1)/2-.002, '\sigma_1', 'fontsize',18)
myplot([0 US(1,2)], [0, US(2,2)],'k--','linewidth',1);
text(US(1,2)/2+.001, US(2,2)/2, '\sigma_2', 'fontsize',18)
xlabel('$x$','Interpreter','Latex')
ylabel('$\dot{x}$','Interpreter','Latex')

axis equal
% 
% I = eye(2);
% 
% vi = V'*I;
% svi = S*vi;
% usvi = U*svi;
% AI = At*I
% % lets look at some some plots
% figure(2), 
% plot([0 I(1,1)], [0 I(2,1)], 'b'), hold on
% plot([0 I(1,2)], [0 I(2,2)], 'b:')
% grid on
% 
% plot([0 vi(1,1)], [0 vi(2,1)], 'r'), hold on
% plot([0 vi(1,2)], [0 vi(2,2)], 'r:'), 
% 
% plot([0 svi(1,1)], [0 svi(2,1)], 'g'), hold on
% plot([0 svi(1,2)], [0 svi(2,2)], 'g:'), 
%  
% plot([0 usvi(1,1)], [0 usvi(2,1)], 'y'), hold on
% plot([0 usvi(1,2)], [0 usvi(2,2)], 'y:'), 
%    
% plot([0 AI(1,1)], [0 AI(2,1)], 'k--'), hold on
% plot([0 AI(1,2)], [0 AI(2,2)], 'k:'), 