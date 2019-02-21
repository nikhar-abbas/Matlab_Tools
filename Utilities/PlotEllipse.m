function PlotEllipse(a,b,theta,vecmat)
% Plots an ellipsoid with major and minor axes of length a and b,
% respectively
%
% Inputs: a - major axis length
%         b - minor axis length
%         theta - rotation of matrix
%         vecmat - addition vectors to plot on ellipse


switch nargin
    case 2
        theta = 0;
        vecmat = [];
    case 3
        vecmat = [];
    case 4
        v1 = vecmat(:,1);
        v2 = vecmat(:,2);
end




t = linspace(0,2*pi,100);

x = a*cos(t)*cos(theta) - b*sin(t)*sin(theta);
y = b*sin(t)*cos(theta) + a*cos(t)*sin(theta);


myplot(x,y,'k');

if ~isempty(vecmat)
    hold on
    myplot([0 v1(1)],[0 v1(2)],'r');
    myplot([0 v2(1)],[0 v2(2)],'r');
end



