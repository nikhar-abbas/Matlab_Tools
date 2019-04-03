function [tvec] = Pre_LinTimes(omega, t0, nlin, nrot)
% Use this function to find the linearization times given a fixed rotor
% speed. 

% This function is for pre-processing use to define the linearization times
% in an openfast input (*.fst) file

% Inputs:
%     omega - rotor speed (rpm)
%     t0 - initial linearization time (s)
%     nlin - number of linearization points per rotation
%     nrot - number of rotations to linearize over
% Outputs:
%     tvec - vector of times for openfast to linearize at
%
% Nikhar Abbas - January 2019
%%
switch nargin
    case 2             % Default to 18 linearization points over 1 rotation
        nlin = 18;
        nrot = 1; 
        display('case = 2')
    case 3             % Default to 1 rotation
        nrot = 1;
        display('case = 1')
end

om_s = omega/60;                    % rot/sec
ts = 1/om_s;                        % sec/rot

tlin = ts/(nlin);                   % sec/nlin


linmat = [0:(nrot*nlin)-1];
tlinmat = linmat.*tlin;
tvec = t0.*ones(1,length(linmat))+tlinmat;      % Output times to linearize

end


