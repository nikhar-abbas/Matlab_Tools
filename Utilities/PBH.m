function [ctr] = PBH(A,B,lam)
% Performs the Popov-Belavitch-Hautus test for controllability
%
% Inputs:   A - A matrix (square)
%           B - B matrix
%           lam - eigenvalues to test
%
% Outputs: ctr - 1 for controllable, 0 for not controllable
%
%

n = length(A);

M = [lam*eye(n)-A  B];

r = rank(M);

% 1 if controllable, 0 if not
if n == r
    ctr = 1;
elseif n ~= r
    ctr = 0;
else
    error('Somehow you broke the PBH test')
end
    
end

% Nikhar Abbas
% Oct 15, 2018


