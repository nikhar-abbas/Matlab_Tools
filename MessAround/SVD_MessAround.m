%% SVD_MessAround
% Just a script to explore the SVD in a simplified, two DOF model.

theta = 45;

% Make a fancy dynamic plot
[Wc, U,S,V, E, Ex] = svd_an(theta)

% Pl_SVD(Wc)
%% Make a fancy dynamic plot
f = gcf;
h = uicontrol('Parent',f,'Style','slider','Position',[400,0,400,40],...
              'value',theta, 'min',0, 'max',90);

h.Callback = @(es,ed) svd_an(es.Value);



%%
function [Wc, U,S,V, E, Ex] = svd_an(theta)
%% State Space
m = 10;
% theta = 45;
theta
k1 = 4;
k2 = k1;
b1 = 4;
b2 = b1;
A = [0 1 0 0;...
    -k1/m -b1/m 0 0;...
    0 0 0 1;...
    0 0 -k2/m -b2/m];
B = [0 0;...
    1/m sind(theta)/m;...
    0 0;...
    0 cosd(theta)/m];
C = eye(4);
D = zeros(4,2);

sys = ss(A,B,C,D);


%% Controllability Gramian, and SVD
Wc = gram(sys,'c');
[U,S,V] = svd(Wc)

%% Controllable Directions and Energy
% Directions
v1 = V(:,1);
v2 = V(:,3);

% SVD Energy
E = V'/Wc*V;
Esvd = diag(E);
% v1E = v1.*Esvd(1);
% v2E = v2.*Esvd(2);
v1E = v1'/Wc*v1
v2E = v2'/Wc*v2
% Singular State Energy
xm = eye(4);
Ex = xm'/Wc*xm;
Ess = diag(Ex)




%% Some plots
f = figure(10);
set(gca,'ColorOrderIndex',1);
myplot([0 v1(1)],[0 v1(3)]); hold on
% myplot([0 v2(1)],[0 v2(2)]);
% myplot([
set(gca,'ColorOrderIndex',1);
% myplot([0 v1E(1)],[0 v1E(3)],'--');
% myplot([0 v2E(1)],[0 v2E(2)],'--');
% 
% myplot(Ess(1),Ess(3),'o'); set(gca,'ColorOrderIndex',3);
% myplot(Ess(1),Ess(3),'o','markersize',12); set(gca,'ColorOrderIndex',3);
% myplot(Ess(1),Ess(3),'o','markersize',18);

title({'Controllable Directions (solid)';...
'Energy to move "one" unit in that direction (dashed)'})

% hold off
end

%% Make plots function
function Pl_SVD(Wc)

%% Controllability Gramian, and SVD
% Wc = gram(sys,'c');
[U,S,V] = svd(Wc);

%% Controllable Directions and Energy
% Directions
v1 = abs(V(:,1));
v2 = abs(V(:,2));

% SVD Energy
E = V'/Wc*V;
Esvd = diag(E);
v1E = v1.*Esvd(1);
v2E = v2.*Esvd(2);

% Singular State Energy
xm = eye(4);
E = xm'/Wc*xm;
Ess = diag(E);




%% Some plots
f = figure(10);
set(gca,'ColorOrderIndex',1);
myplot([0 v1(1)],[0 v1(3)],'r'); hold on
% myplot([0 v2(1)],[0 v2(2)]);
% myplot([
% set(gca,'ColorOrderIndex',1);
% myplot([0 v1E(1)],[0 v1E(3)],'--');
% myplot([0 v2E(1)],[0 v2E(2)],'--');
% 
% myplot(Ess(1),Ess(3),'o'); set(gca,'ColorOrderIndex',3);
% myplot(Ess(1),Ess(3),'o','markersize',12); set(gca,'ColorOrderIndex',3);
% myplot(Ess(1),Ess(3),'o','markersize',18);

title({'Controllable Directions (solid)';...
'Energy to move "one" unit in that direction (dashed)'})
end