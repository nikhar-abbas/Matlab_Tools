%% Post_GetLinear
% - Load and average linearization data from OpenFast
%
% Inputs: Lindir - output directory of linearization cases
%         LinfileBase - base name for OpenFast runs, i.e.
%                       <LinfileBase>.x.lin
%         nlin - number of linearizations. These should generally span an
%                entire rotation of the rotor
% Outputs: linout - output of linearization data
%          linavg - averaged linearization data over all linearization points
%          

% Note: There is a number of extra "features" to this code that may be used
% for linearization analysis, but will change the output data from openfast
% appropriately
% 
% 
% 
% Nikhar Abbas, January 2019

%%
function [linout, linavg] = Post_GetLinear(Lindir, LinfileBase, nlin)

%% Setup paths
addpath(Lindir)

%% Load Files
% clear linout
for ind = 1:nlin
    % Load linearization
    linfile = [Lindir filesep LinfileBase '.' num2str(ind) '.lin'];
%     linout(ind) = ReadFASTLinear(linfile);
%     FileNames{ind} = [OutfileBase '.' num2str(ind) '.lin'];
    FileNames{ind} = linfile;    
end


%% Multiblade Coordinate Transform
% This takes the mixed rotating and non-rotating reference frame outputs
% from the FAST .lin files and transforms them into an entirely nonrotating
% reference frame state space. The idea here "normalizing" the outputs to
% be consistent, and thus enable linear analysis and control design.
%
% All nonrotating frame state space matrices are found using the methods
% that result in equations (29) through (33) of NREL_MBC.pdf (NREL report
% on MBC

% for ind = 1:nlin
%     FileNames = FileNameVec(ind);
%     GetMats_f8;
%     mbc3;
%     
%     if exist('MBC_A')
%         linout(ind).A_MBC = MBC_A;
%         linout(ind).B_MBC = MBC_B;
%         linout(ind).C_MBC = MBC_C;
%         linout(ind).D_MBC = MBC_D;
%     else
%         linout(ind).A_MBC = linout(ind).A;
%         linout(ind).B_MBC = linout(ind).B;
%         linout(ind).C_MBC = linout(ind).C;
%         linout(ind).D_MBC = linout(ind).D;
%     end
%     clear MBC_A MBC_B MBC_C MBC_D
% 
% end


[matData, linout] = fx_getMats(FileNames)
try
    [MBC, matData, FAST_linData] = fx_mbc3(FileNames) 
catch
    display('There are no rotating states - but this is okay, MBC3 is just kind of silly to do.')
end

for ind = 1:nlin
    if exist('MBC')
        linout(ind).A_MBC = MBC.A(:,:,nlin);
        linout(ind).B_MBC = MBC.B(:,:,nlin);
        linout(ind).C_MBC = MBC.C(:,:,nlin);
        linout(ind).D_MBC = MBC.D(:,:,nlin);
    else
        linout(ind).A_MBC = linout(ind).A;
        linout(ind).B_MBC = linout(ind).B(:,1:9);
        linout(ind).C_MBC = linout(ind).C;
        linout(ind).D_MBC = linout(ind).D(:,1:9);
    end
end
    


%% Some analysis


for ind = 1:nlin

% Some organization
    A_MBC = linout(ind).A_MBC;
    B_MBC = linout(ind).B_MBC;
    C_MBC = linout(ind).C_MBC;
    D_MBC = linout(ind).D_MBC;

%% Change some states and input units if desired                -- nja: commented out initially
% % Change units on some specific states 
%     M_x = ones(1,length(linout(ind).x_desc));
%     % change rad to deg
%     rad_st = contains(linout(ind).x_desc,'rad');
%     rad_chstates = find(rad_st == 1)';
%     M_x(rad_chstates) = 180/pi;     
%     linout(ind).x_desc(rad_st) = strrep(linout(ind).x_desc(rad_st),'rad','deg');
%     % change Nm to kNm
%     Nm_st = contains(linout(ind).x_desc,'Nm');
%     Nm_chstates = find(Nm_st == 1)';
%     M_x(Nm_chstates) = 1000;
%     linout(ind).x_desc(rad_st) = strrep(linout(ind).x_desc(rad_st),'Nm','kNm');
%     M_x = diag(M_x);
%     
% 
% % Change units on some specific inputs
%     M_u = ones(1,length(linout(ind).u_desc));
%     % change rad to deg
%     rad_u = contains(linout(ind).u_desc,'rad');
%     rad_chu = find(rad_u == 1)';
%     M_u(rad_chu) = 180/pi;      
%     linout(ind).u_desc(rad_u) = strrep(linout(ind).u_desc(rad_u),'rad','deg');
%     % change Nm to kNm
%     Nm_u = contains(linout(ind).u_desc,'Nm');
%     Nm_chu = find(Nm_u == 1)';
%     M_u(Nm_chu) = 1000;
%     linout(ind).u_desc(Nm_u) = strrep(linout(ind).u_desc(Nm_u),'Nm','kNm');
%     M_u = diag(M_u);
% 
% Apply transformation matrices
%     A_MBC = M_x*A_MBC/M_x;
%     B_MBC = M_x*B_MBC/M_u;
%     C_MBC = C_MBC/M_x;
%     D_MBC = D_MBC/M_u; 
    

% Remove unnecessary/unstable states and inputs                    -- nja: this might be desired for some post-processing/analysis
    gen_st = contains(linout(ind).x_desc,'Variable speed generator DOF');
    genstates = find(gen_st == 1)';                                % remove ED Variable speed generator DOF
    if isempty(genstates) 
        rmstates = [];
    else
        rmstates = [genstates(1)];    
    end
    rminputs = [];                                         % 1:7 = All inputs besides GenTq and ColBldPitch
    A_MBC(rmstates,:) = []; A_MBC(:,rmstates) = [];
    B_MBC(rmstates,:) = []; 
    B_MBC(:,rminputs) = [];
    C_MBC(:,rmstates) = [];
    D_MBC(:,rminputs) = [];
    sys_rm = ss(A_MBC, B_MBC, C_MBC, D_MBC);
    linout(ind).x_desc_rm = linout(ind).x_desc;
    linout(ind).x_desc_rm(rmstates) = [];
    linout(ind).u_desc_rm = linout(ind).u_desc;
    linout(ind).u_desc_rm(rminputs) = [];
   
    
%% Save states matrices
if ~exist('rmstates'), rmstates = []; end
if ~exist('rminputs'), rminputs = []; end

    linout(ind).A_MBC = A_MBC;
    linout(ind).B_MBC = B_MBC;
    linout(ind).C_MBC = C_MBC;
    linout(ind).D_MBC = D_MBC;
    linout(ind).x_desc = linout(ind).x_desc;
    linout(ind).x_desc(rmstates) = [];
    linout(ind).u_desc = linout(ind).u_desc;
    linout(ind).u_desc(rminputs) = [];

    % Eigenvalues of system
    linout(ind).eigs = eig(A_MBC);
    
    % Clear vars
    clear A_MBC B_MBC C_MBC D_MBC M_x M_u 

end

%% Some averaging over linearized states

Asum = [zeros(size(linout(1).A_MBC))];
Bsum = [zeros(size(linout(1).B_MBC))];
Csum = [zeros(size(linout(1).C_MBC))];
Dsum = [zeros(size(linout(1).D_MBC))];
for ind = 1:nlin
   Asum = Asum + linout(ind).A_MBC;
   Bsum = Bsum + linout(ind).B_MBC;
   Csum = Csum + linout(ind).C_MBC;
   Dsum = Dsum + linout(ind).D_MBC;
end

linavg.A = Asum./nlin;
linavg.B = Bsum./nlin;
linavg.C = Csum./nlin;
linavg.D = Dsum./nlin;

linavg.x_desc = linout(1).x_desc;
linavg.u_desc = linout(1).u_desc;
linavg.y_desc = linout(1).y_desc;
%% shorten states descriptions
% st = linavg.x_desc;
% for j = 1:length(st)
%     if contains(st(j),'blade 1')
%         st{j} = strrep(st{j},'blade 1', 'coning mode');
%     elseif contains(st(j),'blade 2')
%         st{j} = strrep(st{j},'blade 2', 'cosine-cyclic mode');
%     elseif contains(st(j),'blade 3')
%         st{j} = strrep(st{j},'blade 3', 'sine-cyclic mode');
%     end
%     if contains(st(j), 'First time derivative')
%         st{j} = strrep(st{j},'First time derivative', 'd/dt');
%     end
%     if contains(st(j), 'DOF')
%         st{j} = strrep(st{j},' DOF', '');
%     end
%     if contains(st(j), 'ED')
%         st{j} = strrep(st{j},'        ED', '');
%     end
%     
%     sta = strfind(st{j},'(');
%     if contains(st(j), 'm/s')
%         en = length(st{j}) - 5;
%     elseif contains(st(j), ["rad/s","deg/s"])
%         en = length(st{j}) - 7;
%     elseif contains(st(j), ["rad","deg"])
%         en = length(st{j}) - 5;
%     elseif contains(st(j), 'm')
%         en = length(st{j}) - 3;
%     end
%    
%     st{j}(sta-1:en) = [];
%     st{j}(1:3) = [];
% end
% 
% linavg.st_desc = st;

end




