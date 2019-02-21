function tau_a = AeroTorque(Cp,lambda,v)
% This function calculates the aerodynamic torque of a defined wind turbine
% rotor
% 
% Inputs: Cp - Coefficient of Power
%         lambda - Tip Speed Ratio
%         v - wind speed (m/s)
%         
% Outputs: tau_a - Aerodynamic Torque (N-m)
%
% Nikhar Abbas, January 2019


%% -- Parameters for DTU 10MW --
J = 156348032.208;                  % Rotor Inertia (kg-m^2)
rho = 1.225;                        % Air Density (kg/m^3)
R = 89.166;                         % Rotor Radius (m)
A = pi*R^2;                         % Rotor Swept Area (m^2)


%% Aerodynamic Torque Calculation 
tau_a = 1/2 * rho*A*R*v^2 * Cp/lambda;




end