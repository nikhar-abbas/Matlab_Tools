%% SVD_MessAround
% Just a script to explore the SVD in a simplified, two DOF model.

a1 = 0;

% svd_an(a1)
% f = gcf;
% h = uicontrol('Parent',f,'Style','slider','Position',[400,0,400,40],...
%               'value',a1, 'min',0, 'max',1);

h.Callback = @(es,ed) svdanal(es.Value);

%%
% function svd_an(a1)
%% State Space
% a2 = 0;
% a2 = get(hObject,'value')
A = [-1 a1; ...
     0 -1]
B = [1 0; ...
     0 1]
C = [1 0; 0 1];
D = zeros(2,2);

sys = ss(A,B,C,D);


%% Controllability Gramian, and SVD

Wc = gram(sys,'c');
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
xm = eye(2);
E = xm'/Wc*xm;
Ess = diag(E);




%% Some plots
f = figure(10);
set(gca,'ColorOrderIndex',1);
myplot([0 v1(1)],[0 v1(2)]); hold on
myplot([0 v2(1)],[0 v2(2)]);
% myplot([
set(gca,'ColorOrderIndex',1);
myplot([0 v1E(1)],[0 v1E(2)],'--');
myplot([0 v2E(1)],[0 v2E(2)],'--');

myplot(Ess(1),Ess(2),'o'); set(gca,'ColorOrderIndex',3);
myplot(Ess(1),Ess(2),'o','markersize',12); set(gca,'ColorOrderIndex',3);
myplot(Ess(1),Ess(2),'o','markersize',18);

title({'Controllable Directions (solid)';...
'Energy to move "one" unit in that direction (dashed)'})

hold off
end