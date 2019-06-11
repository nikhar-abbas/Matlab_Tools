% GenStepWindFile
% Function for generating a step wind speed file given:
% 1) A base filename
% 2) An array of windspeeds
% - The generated file will have a name BaseName_Windspeed_MPS.wnd
%
%In:    fileName          -   filename for step file
%       tstep             -   length of timestep between wind changes
%       windSpeedArray    -   array of wind speeds in m/s
%
%
% Paul Fleming, JUNE 2011
% Using code copied from functions written by Jan-Willem van Wingerden in
% July 2010
%
% May 2019 - Modified by Nikhar Abbas for some added flexibility. Could be
% run using Pre_GenStepWind.m for even more ease!


function GenStepWindFile(fileName,tstep, windSpeed)


fid= fopen(fileName,'w');
fprintf(fid,'!Wind file with step changes in wind speed. \n');
fprintf(fid,'!Time  Wind     Wind	Vert.       Horiz.      Vert.       LinV        Gust \n');
fprintf(fid,'!      Speed    Dir    Speed       Shear		Shear       Shear       Speed \n');
fprintf(fid,'%2.2f %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f\n', 0, windSpeed(1), 0, 0,0, 0,0,0);
for i=1:1:length(windSpeed)
    fprintf(fid,'%2.2f %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f\n', 0.1+tstep*(i-1), windSpeed(i), 0, 0,0, 0,0,0);
    fprintf(fid,'%2.2f %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f\n', tstep*i, windSpeed(i), 0, 0,0, 0,0,0);    
end

fclose(fid);