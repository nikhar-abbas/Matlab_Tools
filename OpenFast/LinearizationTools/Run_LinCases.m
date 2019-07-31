%% Run_LinCases

% This script is written to cycle through a number of cases for
% linearized model analysis. Degrees of freedom are switched on or off in
% the ElastDyn module, and a corresponding steady state BldPitch operating
% point is prescribed. 
%
% NOTE: On OSX, this script needs to be run via the terminal - or via a
% MATLAB gui initiated in the terminal - due to $PATH configurations.
%
%
%
% Setup: The LinTimes input in the openfast input file should be
% initialized, and the 'Linearize' input should be switched to 'True'.
% A corresponding case file for each case defined in this script should be 
% made. This case file should include any/all files necessary to run a
% complete openfast simulation with linearization on.  This file should be
% titled 'CaseX', where 'X' is the case number. The openfast files should
% have the format '(modelname)_(module).dat' and '(modelname).fst'. I.E.
% 'DTU_10MW.fst', with 'DTU_10MW_ElastoDyn.dat'. 

% Running: Define all degrees of freedom available in openfast (listed
%        below). Then define non-floating DOF options. If floating is
%        enabled, the desired floating degrees of freedom will be appended
%        to the fixed bottom cases. 
%        *****
%        Make sure cpath is properly defined to point to the right
%        (numbered) case folders 
%        *****

% DOF Options: (at the time of writing this script)
% 	FlapDOF1 
% 	FlapDOF2 
% 	EdgeDOF  
% 	TeetDOF  
% 	DrTrDOF  
% 	GenDOF   
% 	YawDOF   
% 	TwFADOF1 
% 	TwFADOF2 
% 	TwSSDOF1 
% 	TwSSDOF2 
% 	PtfmSgDOF
% 	PtfmSwDOF
% 	PtfmHvDOF
% 	PtfmRDOF 
% 	PtfmPDOF 


% Nikhar Abbas - February 2019

%% Some Initialization Housekeeping
LinFolder = '/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU10MWRWT_NAUTILUS_GoM_FAST_v1.00/Linearizations/';
BaseFile = 'DTU_10MW_NAUTILUS_GoM';

% Number of cases
cnum = 8;       

% Floating?
floating = 1;

%% Define DOF options:
DOFopt = {'FlapDOF1',... 
	'FlapDOF2',... 
	'EdgeDOF',...  
	'TeetDOF',...  
	'DrTrDOF',...  
	'GenDOF',...   
	'YawDOF',...   
	'TwFADOF1',... 
	'TwFADOF2',... 
	'TwSSDOF1',... 
	'TwSSDOF2',... 
	'PtfmSgDOF',...
	'PtfmSwDOF',...
	'PtfmHvDOF',...
	'PtfmRDOF',... 
	'PtfmPDOF',...
    'PtfmYDOF'};

%% Initialize cases

% Case 1
cases.c1.wind = 14;
cases.c1.dof = {'GenDOF'};
% Case 2
cases.c2.wind = 14;
cases.c2.dof = {'GenDOF','TwFADOF1'};
% Case3 
cases.c3.wind = 14;
cases.c3.dof = {'GenDOF', 'TwFADOF1', 'YawDOF'}; 	
% Case 4
cases.c4.wind = 14;
cases.c4.dof = {'GenDOF', 'TwFADOF1', 'YawDOF', 'DrTrDOF'};	 
% Case 5
cases.c5.wind = 14;
cases.c5.dof = {'GenDOF', 'TwFADOF1', 'YawDOF', 'DrTrDOF', 'TwSSDOF1'};		 
% Case 6
cases.c6.wind = 14;
cases.c6.dof = {'GenDOF', 'TwFADOF1', 'YawDOF', 'DrTrDOF', 'TwSSDOF1', 'FlapDOF1'};
% Case 7
cases.c7.wind = 14;
cases.c7.dof = {'GenDOF', 'TwFADOF1', 'YawDOF', 'DrTrDOF', 'TwSSDOF1', 'FlapDOF1', 'EdgeDOF'};
% Case 8
cases.c8.wind = 14;
cases.c8.dof = {'GenDOF', 'TwFADOF1', 'YawDOF', 'DrTrDOF', 'TwSSDOF1', 'FlapDOF1', 'EdgeDOF', 'FlapDOF1', 'FlapDOF2', 'TwFADOF2', 'TwSSDOF2'};

% Add Floating Cases
fcases = {'PtfmPDOF','PtfmSgDOF'};
if floating
    for ci = 1:cnum
        cases.(['c',num2str(ci)]).dof = [cases.(['c',num2str(ci)]).dof, fcases];
    end
end
            


%% Run Setup and Linearization Simulations

for ci = 1:cnum
    
    % Setup file names paths
    cpath = ['Case', num2str(ci),'_PitSur'];
    elastofile = [LinFolder filesep cpath filesep BaseFile '_ElastoDyn.dat'];
    servofile = [LinFolder filesep cpath filesep BaseFile '_ServoDyn.dat'];
    modfile = [LinFolder filesep cpath filesep BaseFile '.fst'];
    outfile = [LinFolder filesep cpath filesep BaseFile '.out'];
    
    % Modify (model).fst
    dt = 0.025;
    TMax = 550;
    lin = 'False';
    mparams = {'DT', 'TMax', 'Linearize'};
    mvals = {dt, TMax, lin};
    Pre_FastMod(modfile,mparams,mvals);
%     
    % Modify Elastodyn
    DOF_flags = cell(length(DOFopt),1);
    DOF_flags(:) = {'False'};
    [~,ia,~] = intersect(DOFopt, cases.(['c', num2str(ci)]).dof);
    DOF_flags(ia) = {'True'};
    Pre_FastMod(elastofile, DOFopt, DOF_flags);    
%     
%     % Modify ServoDyn         ----------------         % Parts of this can be commented out if ElastoDyn's BldPitch(x) is initialized properly alreadly  
    cparams = {'PCMode', 'VSContrl'};
    cvals = {5, 5};
    Pre_FastMod(servofile, cparams, cvals);
    
    % Run Initialization Simulation
    cd([LinFolder filesep cpath])
    system(['openfast ' modfile])
    fo_init = Post_LoadFastOut(outfile);
    
    % Re-define input files                                             
    ss_bldpitch = mean(fo_init.BldPitch1(floor(end*(3/4)),end));            % average BldPitch over last 3/4 of init simulation
    Pre_FastMod(modfile,{'Linearize'},{'True'})
    Pre_FastMod(servofile,cparams,{0,1});
    Pre_FastMod(elastofile,{'BlPitch'},{ss_bldpitch});
    
    % Run Linearization
    system(['openfast ' modfile])
    
end
        
        
        
        
        
%% Find Op Point
sig = 'PtfmYaw'

op = (max(simout.(sig)(3*end/4:end))+min(simout.(sig)(3*end/4:end)))/2

        
        
        
        
        
        
        
        
        
        