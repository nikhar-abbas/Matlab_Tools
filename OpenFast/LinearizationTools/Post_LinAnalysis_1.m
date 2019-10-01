%% Post_LinAnalysis_1
% - Load and analyse a number of FAST linearization files
%
% Notes: Smart uncommenting can let this script stand alone. Good luck...
% 
% Inputs: Outdir - output directory of linearization cases
%         OutfileBase - base name for OpenFast runs
%         nlin - number of linearizations
% Outputs: linout - output of linearization data
%          linavg - averaged linearization data over all linearization points
%          


function [linout, linavg] = Post_LinAnalysis_1(Outdir, OutfileBase, nlin)

%% Setup Paths

addpath(Outdir)

%% Define FileNames
for ind = 1:nlin
    linfile = [Outdir filesep OutfileBase '.' num2str(ind) '.lin'];
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
    
% Find uncrontrollable modes
    lambdas = eig(A_MBC);
    lamrec = [];
    for k = 1:length(lambdas)
        ctr = PBH(A_MBC,B_MBC,lambdas(k));
        lamrec(k) = ctr;
    end
    linout(ind).lambdas = lambdas;
    linout(ind).lamPBHres = lamrec;             % 0 if uncontrollable
    
% Change units on some specific states 
    M_x = ones(1,length(linout(ind).x_desc));
    % change rad to deg
    rad_st = contains(linout(ind).x_desc,'rad');
    rad_chstates = find(rad_st == 1)';                                 % states in radians
    M_x(rad_chstates) = 180/pi;     
    linout(ind).x_desc(rad_st) = strrep(linout(ind).x_desc(rad_st),'rad','deg');
    % change Nm to kNm
    Nm_st = contains(linout(ind).x_desc,'Nm');
    Nm_chstates = find(Nm_st == 1)';
    M_x(Nm_chstates) = 1000;
    linout(ind).x_desc(rad_st) = strrep(linout(ind).x_desc(rad_st),'Nm','kNm');
%     % normalize by tower height
%     twrheight = 115.63;
%     twr_st = contains(linout(ind).x_desc,', m') + contains(linout(ind).x_desc,'tower');
%     twr_chstates = find(twr_st == 2)';
%     M_x(twr_chstates) = M_x(twr_chstates)/twrheight;
%     % normalize by blade length 
%     bldlength = 86.4;
%     bld_st = contains(linout(ind).x_desc,', m') + contains(linout(ind).x_desc,'blade');
%     bld_chstates = find(bld_st == 2)';
%     M_x(bld_chstates) = M_x(bld_chstates)/bldlength;

    M_x = diag(M_x);
    

% Change units on some specific inputs
    M_u = ones(1,length(linout(ind).u_desc));
    % change rad to deg
    rad_u = contains(linout(ind).u_desc,'rad');
    rad_chu = find(rad_u == 1)';                                 % states in radians
    M_u(rad_chu) = 180/pi;      
    linout(ind).u_desc(rad_u) = strrep(linout(ind).u_desc(rad_u),'rad','deg');
    % change Nm to kNm
    Nm_u = contains(linout(ind).u_desc,'Nm');
    Nm_chu = find(Nm_u == 1)';
    M_u(Nm_chu) = 1000;
    linout(ind).u_desc(Nm_u) = strrep(linout(ind).u_desc(Nm_u),'Nm','kNm');
    % Remove HydroDyn Input, if it exists
    Hd_u = contains(linout(ind).u_desc,'HD');
    Hd_ru = find(Hd_u == 1)';
    M_u(Hd_ru) = [];
    % normalize by tower height
%     twrheight = 115.63;
%     twr_u = contains(linout(ind).u_desc,', m') + contains(linout(ind).u_desc,'tower');
%     twr_chu = find(twr_u == 2)';
%     M_u(twr_chu) = M_u(twr_chu)/twrheight;
%      % normalize by blade length
%     bldlength = 86.4;
%     bld_u = contains(linout(ind).u_desc,', m') + contains(linout(ind).u_desc,'tower');
%     bld_chu = find(bld_u == 2)';
%     M_u(bld_chu) = M_u(bld_chu)/bldlength;   

    % Scale platform pitching
    twrheight = 115.63;
    twr_st = contains(linout(ind).x_desc,'Platform');
    twr_chstates = find(twr_st == 1)';
    M_x(twr_chstates) = M_x(twr_chstates)*10;    

    M_u = diag(M_u);

% Apply transformation matrices
    A_MBC = M_x*A_MBC/M_x;
    B_MBC = M_x*B_MBC/M_u;
    C_MBC = C_MBC/M_x;
    D_MBC = D_MBC/M_u;

     
%     twrheight = 1; 115.63;
%     m_st_T = contains(linout(ind).x_desc,', m') + contains(linout(ind).x_desc,'tower');
%     m_chstates_T = find(m_st_T > 1)';                                      % tower height states
%     A_MBC(m_chstates_T,:) = A_MBC(m_chstates_T,:).* 1/twrheight;             % 1/twr height
%     A_MBC(:,m_chstates_T) = A_MBC(:,m_chstates_T).* twrheight;             % times twr height to keep system behavior
%     B_MBC(m_chstates_T,:) = B_MBC(m_chstates_T,:).* twrheight;
%     
%     bllength = 1; 86.4;
%     m_st_bl = contains(linout(ind).x_desc,', m') + contains(linout(ind).x_desc,'blade');
%     m_chstates_bl = find(m_st_T > 1)';                                      % bllength height states
%     A_MBC(m_chstates_bl,:) = A_MBC(m_chstates_bl,:).* 1/bllength;             % 1/bllength
%     A_MBC(:,m_chstates_bl) = A_MBC(:,m_chstates_bl).* bllength;             % times bllength to keep system behavior
%     B_MBC(m_chstates_bl,:) = B_MBC(m_chstates_bl,:).* bllength;    
    
    % Remove uncontrollable states and unnecessary inputs
    gen_st = contains(linout(ind).x_desc,'ED Variable speed generator DOF');
    gen_st = gen_st + contains(linout(ind).x_desc,'Platform yaw rotation DOF');
    gen_st = gen_st + contains(linout(ind).x_desc,'translation DOF');    
%     gen_st = gen_st + contains(linout(ind).x_desc,'blade');
    genstates = find(gen_st == 1)';                                % remove ED Variable speed generator DOF
    if isempty(genstates) 
        rmstates = [];
    else
        rmstates = [genstates];    
    end
%     rminputs = [1:7];                                         % 1-7 = All inputs besides GenTq and ColBldPitch
    rminputs =[1:3, 7,9];
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
     
    
    % Save removed states
    linout(ind).A_MBC_rm = A_MBC;
    linout(ind).B_MBC_rm = B_MBC;
    linout(ind).C_MBC_rm = C_MBC;
    linout(ind).D_MBC_rm = D_MBC;
    % Try to calculate the controllability Gramian and SVD
    try
        % Gramian
        linout(ind).Gc = gram(sys_rm,'c');
        % Min energy
%         linout(ind).Emin = inv(linout(ind.Gc));
        % SVD of min energy matrix
        [U,S,V] = svd(linout(ind).Gc);
        linout(ind).SVD.U = U;
        linout(ind).SVD.S = diag(S);
        linout(ind).SVD.V = V;
        
        
    catch
        display(['Controllability Gramian could not be calculated for linearization point ', num2str(ind)]);
        
        % Fill the structure if its not controllable
        linout(ind).Gc = [];
        linout(ind).SVD.U = [];
        linout(ind).SVD.S = [];
        linout(ind).SVD.V = [];
        
        % Processing
        linout(ind).ctrldir = [];
        linout(ind).SVD.Xwork = [];
    end
    
    % Eigenvalues
    linout(ind).eigs = eig(A_MBC);
    
    % Clear vars
%     clear A_MBC B_MBC C_MBC D_MBC M_x M_u 

end

% contd = [];
% for j= 1:length(linout)
%     contd = [contd linout(j).contd']
% end
%% Some averaging analysis - CLEAN THIS UP AND SEPARATE IT INTO OTHER FUNCTIONS

% Asum = [zeros(size(linout(1).A_MBC_rm))];
% Bsum = [zeros(size(linout(1).B_MBC_rm))];
% Csum = [zeros(size(linout(1).C_MBC_rm))];
% Dsum = [zeros(size(linout(1).D_MBC_rm))];
% for ind = 1:nlin
%    Asum = Asum + linout(ind).A_MBC_rm;
%    Bsum = Bsum + linout(ind).B_MBC_rm;
%    Csum = Csum + linout(ind).C_MBC_rm;
%    Dsum = Dsum + linout(ind).D_MBC_rm;
% end
% 
% linavg.A = Asum./nlin;
% linavg.B = Bsum./nlin;
% linavg.C = Csum./nlin;
% linavg.D = Dsum./nlin;
linavg.A = linout.A_MBC_rm;
linavg.B = linout.B_MBC_rm;
linavg.C = linout.C_MBC_rm;
linavg.D = linout.D_MBC_rm;

sys_avg = ss(linavg.A,linavg.B,linavg.C,linavg.D);
linavg.x_desc = linout(1).x_desc_rm;

% Gramian
linavg.Gc = gram(sys_avg,'c');
 
% SVD
[U,S,V] = svd(linavg.Gc);
linavg.SVD.U = U;
linavg.SVD.S = diag(S);
linavg.SVD.V = V;

% SVD - Energy
[U,S,V] = svd(inv(linavg.Gc));
linavg.SVDe.U = U;
linavg.SVDe.S = diag(S);
linavg.SVDe.V = V;



%% shorten states descriptions
st = linavg.x_desc;
for j = 1:length(st)
    if contains(st(j),'blade 1')
        st{j} = strrep(st{j},'blade 1', 'coning mode');
    elseif contains(st(j),'blade 2')
        st{j} = strrep(st{j},'blade 2', 'cosine-cyclic mode');
    elseif contains(st(j),'blade 3')
        st{j} = strrep(st{j},'blade 3', 'sine-cyclic mode');
    end
    if contains(st(j), 'First time derivative')
        st{j} = strrep(st{j},'First time derivative', 'd/dt');
    end
    if contains(st(j), 'DOF')
        st{j} = strrep(st{j},' DOF', '');
    end
    if contains(st(j), 'ED')
        st{j} = strrep(st{j},'        ED', '');
    end
    
    sta = strfind(st{j},'(');
    if contains(st(j), 'm/s')
        en = length(st{j}) - 5;
    elseif contains(st(j), ["rad/s","deg/s"])
        en = length(st{j}) - 7;
    elseif contains(st(j), ["rad","deg"])
        en = length(st{j}) - 5;
    elseif contains(st(j), 'm')
        en = length(st{j}) - 3;
    end

    
    
    st{j}(sta-1:en) = [];
    st{j}(1:3) = [];
%     st{j} = strcat(st{j},' - ',num2str(j));
end

linavg.st_desc = st;

end




