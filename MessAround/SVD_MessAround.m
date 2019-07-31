%% SVD_MessAround
% Just a script to explore the SVD in a simplified, two DOF model.

theta = 0;
pl = 2;  % plot?
t = .001;
% Do some analysis, make a plot?
[Wc, U,S,V, A, E, Ex] = svd_an(theta, t, pl)

% Pl_SVD(Wc)
%% Make a fancy dynamic plot
f = gcf;
h = uicontrol('Parent',f,'Style','slider','Position',[400,0,400,40],...
              'value',theta, 'min',0, 'max',90);

h.Callback = @(es,ed) svd_an(es.Value, t, pl);

% %% Make a surface of Energy to get "Home"
% theta = 45;
% pl = 0;
% t = 10;                     % Time to return (s)
% [Wc, U,S,V, A, E, Ex] = svd_an(theta, t, pl);
% 
% gsize = 51;
% sfc = ones(gsize,gsize);          % Make surface of desired location (want to be odd, so middle point exists
% [xs,ys] = size(sfc);   
% centi = [xs/2+1 ys/2+1];    % Center index - this will not be the center if gsize isn't odd
% 
% for xi = 1:xs
%     for yi = 1:ys
%         xloc = (xi-centi(1))/xs;
%         yloc = (yi-centi(2))/ys;
%         
%         x0 = [xloc 0 yloc 0]';           % Starting location
%         xdes = [0 0 0 0]';                 % Final Location
%         
%         E = [expm(A*t)*x0 - xdes]'/Wc*[expm(A*t)*x0 - xdes];
%         
%         Esurf(xi,yi) = E;
%     end
% end
% 
% 
% surf(Esurf)
% view(0,90)


%%
function [Wc, U,S,V, A, E, Ex] = svd_an(theta,t, pl)
%% State Space
m = 10;
% theta = 45;
theta
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


%% Controllability Gramian, and SVD
opt = gramOptions('TimeIntervals',[0 t]);
Wc = gram(sys,'c',opt);
% Wc = gram(sys,'c');
[U,S,V] = svd(Wc);

%% Controllable Directions and Energy
% Directions
vi = find(abs(V(1,:)) > 0);
if length(vi) > 1
    v2i = vi(2);
else
    v2i = vi(1);
end
v1 = V(:,1);
v2 = V(:,v2i);

% SVD Energy
E = V'/Wc*V;
Esvd = diag(E);
v1E = v1.*Esvd(1);
v2E = v2.*Esvd(2);
% v1E = v1'/Wc*v1
% v2E = v2'/Wc*v2
% Singular State Energy
xm = eye(4);
Ex = xm'/Wc*xm;
Ess = diag(Ex);

% Singular Values
sigma = diag(S);


%% Some plots
if pl == 1
    f = figure(10);
    % subplot(2,1,1)
    set(gca,'ColorOrderIndex',1);
    myplot([0 v1(1)*sigma(1)],[0 v1(3)*sigma(1)]); hold on
    myplot([0 v2(1)*sigma(v2i)],[0 v2(3)*sigma(v2i)]);
    % subplot(2,1,2)
    % % set(gca,'ColorOrderIndex',1);
    % % myplot([0 v1E(1)],[0 v1E(3)],'--');, hold on
    % % myplot([0 v2E(1)],[0 v2E(2)],'--');
    % % 
    % myplot(Ess(1),Ess(3),'o'); set(gca,'ColorOrderIndex',3);, hold on
    % myplot(Ess(1),Ess(3),'o','markersize',12); set(gca,'ColorOrderIndex',3);
    % myplot(Ess(1),Ess(3),'o','markersize',18);

    title({'Controllable Directions (solid)'});...
    % 'Energy to move "one" unit in that direction (dashed)'})

%     hold off
    
elseif pl == 2
    %% Make a surface of Energy to get "Home"
    gsize = 51;
    sfc = ones(gsize,gsize);          % Make surface of desired location (want to be odd, so middle point exists
    [xs,ys] = size(sfc);   
    centi = [xs/2+1 ys/2+1];    % Center index - this will not be the center if gsize isn't odd

    for xi = 1:xs
        for yi = 1:ys
            xloc = (xi-centi(1))/xs;
            yloc = (yi-centi(2))/ys;
            xv(xi) = xloc;
            yv(yi) = yloc;
            x0 = [xloc 0 yloc 0]';           % Starting location
            xdes = [0 0 0 0]';                 % Final Location

            E = [expm(A.*t)*x0 - xdes]'/Wc*[expm(A.*t)*x0 - xdes];

            Esurf(xi,yi) = E;
        end
    end

    figure(11)
    surf(Esurf')
    view(2)
    colorbar
    xlabel('X')
    ylabel('Y')
end

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