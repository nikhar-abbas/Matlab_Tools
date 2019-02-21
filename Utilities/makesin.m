function [x,y] = makesin(A,f,T,psi)
% Make a sine wave

% Inputs: A - Amplitude
%         f = frequency
%         T = Total Time
% Outputs: x = Time vector
%          y = Sine Wave

switch nargin
    case 2
        T = 1;
        psi = 0;
    case 3
        psi = 0;
end 

% A = 1;
% f = .5;
% T = 100;

x = 0:pi/100:T;
y = A*sin(2*pi*f*x + psi);

end