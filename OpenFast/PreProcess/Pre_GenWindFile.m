 % GenStepWindFile
% Function for generating a wind speed file given:
% 1) A base filename
% 2) An array of windspeeds
%
%In:    fileName          -   filename for output file
%       time              -   array of timesteps
%       windspeed         -   array of wind speeds in m/s
%
%
% Paul Fleming, JUNE 2011
% Using code copied from functions written by Jan-Willem van Wingerden in
% July 2010
%
% May 2019 - Modified by Nikhar Abbas for some added flexibility. Could be
% run using Pre_GenStepWind.m for even more ease!


function Pre_GenWindFile(fileName,time, windspeed)


fid= fopen(fileName,'w');
fprintf(fid,'!Wind file for input to OpenFAST. \n');
fprintf(fid,'!Time  Wind     Wind	Vert.       Horiz.      Vert.       LinV        Gust \n');
fprintf(fid,'!      Speed    Dir    Speed       Shear		Shear       Shear       Speed \n');
fprintf(fid,'%2.2f %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f\n', 0, windspeed(1), 0, 0,0, 0,0,0);
for i=1:1:length(windspeed)
    fprintf(fid,'%2.2f %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f\n', time(i), windspeed(i), 0, 0,0, 0,0,0);
    fprintf(fid,'%2.2f %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f\n', time(i), windspeed(i), 0, 0,0, 0,0,0);    
end

fclose(fid);