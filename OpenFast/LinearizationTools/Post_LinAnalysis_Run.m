%% Post_LinAnalysis_Run
% Script to run Post_LinAnalysis_1 and make some fancy plots


%% Plot and Save?
mkplot = 1;
savepl = 0;

savedir = '/Users/nabbas/Documents/Notes/ResearchNotes/Figures/';

%% Load Linearization data

% define paths
% Lindir = '/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU10MWRWT/Linearizations/AboveRated/';
% OutfileBase = 'DTU_10MW_RWT';
Lindir = '/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU10MWRWT_NAUTILUS_GoM_FAST_v1.00/Linearizations/';
OutfileBase = 'DTU_10MW_NAUTILUS_GoM';
nlin = 24;
cnum = 7;

% Load cases
for cind = 1:cnum
    cdir = ['Case' num2str(cind)];
    [linout, linavg] = Post_LinAnalysis_1([Lindir, cdir], OutfileBase, nlin);
    
    lindata(cind).linout = linout;
    lindata(cind).linavg = linavg;
end

[FS_desc, SS_desc] = getstates;

%% Process the data
contmat = zeros(length(FS_desc),cnum);

for cind = 1:cnum
    Vdat = lindata(cind).linavg.SVD.V(:,1);
    cvec = Vdat;
%     cvec = lindata(cind).linavg.SVDe.S(1)*Vdat;

    % Find state location in list - this will break if getstates is not the
    % same as the state name simplification in Post_LinAnalysis_1
    [~,st_loc] = ismember(lindata(cind).linavg.st_desc, SS_desc);
    
    % Energy
    xdes = eye(size(lindata(cind).linavg.Gc));
    Emat = xdes/lindata(cind).linavg.Gc*xdes;
    cvec_e = diag(Emat);
    
    % fill the matrix of controllable direction
    contmat(st_loc,cind) = abs(cvec);
    contmat_E(st_loc,cind) = abs(cvec_e);
%     contmat = eye(size(lindata(cind).lin
end

zrws = find(sum(contmat')' ~= 0);

contmat = contmat(zrws,:);
contmat_E = contmat_E(zrws,:);
SS_desc_f = SS_desc(zrws);


%% Make your fancy plots
if mkplot
    fi=1;           % change this to change initial figure
    fig(1) = figure(fi);
    gh = bar3(abs(contmat));
    view(2)
    for k = 1:length(gh)
        zdata = gh(k).ZData;
        gh(k).CData = zdata;
        gh(k).FaceColor = 'interp';
    end
    
    myColorMap = jet(256);
    myColorMap(1,:) = 1;
    colormap(myColorMap);
    colorbar('southoutside')
    
    ax = gca;
    ax.YTick = [1:length(SS_desc_f)];
    ax.YTickLabel = SS_desc_f;
    ax.CLim = [0 0.1];
    ax.FontSize = 9; 
    
    xlabel('Case Number')
    title('Color Scaled Controllable Directions','FontSize',10,'horizontalAlignment', 'right')
    
    
    fig(2) = figure(fi+1);
    gh = bar3(abs(contmat));
    view(2)
    for k = 1:length(gh)
        zdata = gh(k).ZData;
        gh(k).CData = zdata;
        gh(k).FaceColor = 'interp';
    end
    
    myColorMap = jet(256);
    myColorMap(1,:) = 1;
    colormap(myColorMap);
    colorbar('southoutside')
    
    ax = gca;
    ax.YTick = [1:length(SS_desc_f)];
    ax.YTickLabel = SS_desc_f;
    ax.FontSize = 9; 
    
    xlabel('Case Number')
    title('Controllable Directions','FontSize',10,'horizontalAlignment', 'right')
    
    fig(3) = figure(fi+2);
    Eplot1 = log10(abs(contmat_E));
    ninf = find(Eplot1==-inf);
    Eplot = abs(Eplot1);
    Eplot(ninf) = 0;
    gh = bar3(Eplot);
    view(2)
    for k = 1:length(gh)
        zdata = gh(k).ZData;
        gh(k).CData = zdata;
        gh(k).FaceColor = 'interp';
    end
    
    myColorMap = jet(256);
    myColorMap(1,:) = 1;
    colormap(myColorMap);
    colorbar('southoutside')
    
    ax = gca;
    ax.YTick = [1:length(SS_desc_f)]; 
    ax.YTickLabel = SS_desc_f;
    ax.FontSize = 9; 
    
    xlabel('Case Number')
    title('Control Energy','FontSize',10,'horizontalAlignment', 'right')
    
    if savepl
        s_str = {'ContDir_scl', 'ContDir','ContEng'};
        
        for nfig = 1:length(fig)
            fig(nfig).Units = 'inches';
            fig(nfig).PaperPosition = [0 0 6.5 5];
            saveas(fig(nfig),[savedir s_str{nfig}],'epsc')
        end
    end
end

%% Get state names
% mostly putting this in a function to keep processing code clean 

function [FS_desc, SS_desc] = getstates
FS_desc = {'ED Platform horizontal surge translation DOF (internal DOF index = DOF_Sg), m';...
        'ED Platform horizontal sway translation DOF (internal DOF index = DOF_Sw), m';...
        'ED Platform vertical heave translation DOF (internal DOF index = DOF_Hv), m';...
        'ED Platform roll tilt rotation DOF (internal DOF index = DOF_R), rad';...
        'ED Platform pitch tilt rotation DOF (internal DOF index = DOF_P), rad';...
        'ED Platform yaw rotation DOF (internal DOF index = DOF_Y), rad';...
        'ED 1st tower fore-aft bending mode DOF (internal DOF index = DOF_TFA1), m';...
        'ED 1st tower side-to-side bending mode DOF (internal DOF index = DOF_TSS1), m';...
        'ED 2nd tower fore-aft bending mode DOF (internal DOF index = DOF_TFA2), m';...
        'ED 2nd tower side-to-side bending mode DOF (internal DOF index = DOF_TSS2), m';...
        'ED Nacelle yaw DOF (internal DOF index = DOF_Yaw), rad';...
        'ED Variable speed generator DOF (internal DOF index = DOF_GeAz), rad';...
        'ED Drivetrain rotational-flexibility DOF (internal DOF index = DOF_DrTr), rad';...
        'ED 1st flapwise bending-mode DOF of blade 1 (internal DOF index = DOF_BF(1,1)), m';...
        'ED 1st flapwise bending-mode DOF of blade 2 (internal DOF index = DOF_BF(2,1)), m';...
        'ED 1st flapwise bending-mode DOF of blade 3 (internal DOF index = DOF_BF(3,1)), m';...
        'ED 1st edgewise bending-mode DOF of blade 1 (internal DOF index = DOF_BE(1,1)), m';...
        'ED 1st edgewise bending-mode DOF of blade 2 (internal DOF index = DOF_BE(2,1)), m';...
        'ED 1st edgewise bending-mode DOF of blade 3 (internal DOF index = DOF_BE(3,1)), m';...
        'ED 2nd flapwise bending-mode DOF of blade 1 (internal DOF index = DOF_BF(1,2)), m';...
        'ED 2nd flapwise bending-mode DOF of blade 2 (internal DOF index = DOF_BF(2,2)), m';...
        'ED 2nd flapwise bending-mode DOF of blade 3 (internal DOF index = DOF_BF(3,2)), m';...
        'ED First time derivative of Platform horizontal surge translation DOF (internal DOF index = DOF_Sg), m/s';...
        'ED First time derivative of Platform horizontal sway translation DOF (internal DOF index = DOF_Sw), m/s';...
        'ED First time derivative of Platform vertical heave translation DOF (internal DOF index = DOF_Hv), m/s';...
        'ED First time derivative of Platform roll tilt rotation DOF (internal DOF index = DOF_R), rad/s';...
        'ED First time derivative of Platform pitch tilt rotation DOF (internal DOF index = DOF_P), rad/s';...
        'ED First time derivative of Platform yaw rotation DOF (internal DOF index = DOF_Y), rad/s';...
        'ED First time derivative of 1st tower fore-aft bending mode DOF (internal DOF index = DOF_TFA1), m/s';...
        'ED First time derivative of 1st tower side-to-side bending mode DOF (internal DOF index = DOF_TSS1), m/s';...
        'ED First time derivative of 2nd tower fore-aft bending mode DOF (internal DOF index = DOF_TFA2), m/s';...
        'ED First time derivative of 2nd tower side-to-side bending mode DOF (internal DOF index = DOF_TSS2), m/s';...
        'ED First time derivative of Nacelle yaw DOF (internal DOF index = DOF_Yaw), rad/s';...
        'ED First time derivative of Variable speed generator DOF (internal DOF index = DOF_GeAz), rad/s';...
        'ED First time derivative of Drivetrain rotational-flexibility DOF (internal DOF index = DOF_DrTr), rad/s';...
        'ED First time derivative of 1st flapwise bending-mode DOF of blade 1 (internal DOF index = DOF_BF(1,1)), m/s';...
        'ED First time derivative of 1st flapwise bending-mode DOF of blade 2 (internal DOF index = DOF_BF(2,1)), m/s';...
        'ED First time derivative of 1st flapwise bending-mode DOF of blade 3 (internal DOF index = DOF_BF(3,1)), m/s';...
        'ED First time derivative of 1st edgewise bending-mode DOF of blade 1 (internal DOF index = DOF_BE(1,1)), m/s';...
        'ED First time derivative of 1st edgewise bending-mode DOF of blade 2 (internal DOF index = DOF_BE(2,1)), m/s';...
        'ED First time derivative of 1st edgewise bending-mode DOF of blade 3 (internal DOF index = DOF_BE(3,1)), m/s';...
        'ED First time derivative of 2nd flapwise bending-mode DOF of blade 1 (internal DOF index = DOF_BF(1,2)), m/s';...
        'ED First time derivative of 2nd flapwise bending-mode DOF of blade 2 (internal DOF index = DOF_BF(2,2)), m/s';...
        'ED First time derivative of 2nd flapwise bending-mode DOF of blade 3 (internal DOF index = DOF_BF(3,2)), m/s';...
        };

% shorten states descriptions
% ***** Note - this needs to be the same script as is in Post_LinAnalysis_1,
% otherwise you might break things

sts = FS_desc;
for j = 1:length(sts)
    if contains(sts(j),'blade 1')
        sts{j} = strrep(sts{j},'blade 1', 'coning mode');
    elseif contains(sts(j),'blade 2')
        sts{j} = strrep(sts{j},'blade 2', 'cosine-cyclic mode');
    elseif contains(sts(j),'blade 3')
        sts{j} = strrep(sts{j},'blade 3', 'sine-cyclic mode');
    end
    if contains(sts(j), 'First time derivative')
        sts{j} = strrep(sts{j},'First time derivative', 'd/dt');
    end
    if contains(sts(j), 'DOF')
        sts{j} = strrep(sts{j},' DOF', '');
    end
    if contains(sts(j), 'ED')
        sts{j} = strrep(sts{j},'        ED', '');
    end
    
    sta = strfind(sts{j},'(');
    if contains(sts(j), 'm/s')
        en = length(sts{j}) - 5;
    elseif contains(sts(j), ["rad/s","deg/s"])
        en = length(sts{j}) - 7;
    elseif contains(sts(j), ["rad","deg"])
        en = length(sts{j}) - 5;
    elseif contains(sts(j), 'm')
        en = length(sts{j}) - 3;
    end

    % replace rad with deg
    if contains(sts(j), 'rad')     
        sts(j) = strrep(sts(j),'rad','deg');
    end
    sts{j}(sta-1:en) = [];
    sts{j}(1:3) = [];
%     st{j} = strcat(st{j},' - ',num2str(j));
end



SS_desc = sts;

end










