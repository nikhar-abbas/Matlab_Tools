%% PlotOutputs.m
% Plots desired outputs from defined OpenFAST .out file or simulink s-function output 

%% Load fAST data

% simulink or fast?
type = 1;                   % 1 for fast outfile, 2 for simulink 
% save figures?
sfig = 0;                   % 1 to save figures
%save data?
sdata = 0;                  % 1 to save data

if type == 1
    % Load FAST output file
%     dir = '/Users/nabbas/Google Drive File Stream/My Drive/PaoLabShare/FloatingStudy/TurbineModels/NREL5MW-OC3Spar/00-Base';
%     filename = 'Model_u18.out';
%     file = '/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU10MWRWT/Linearizations/LinTest/DTU_10MW_RWT.out';
    dir = '/Users/nabbas/Google Drive File Stream/My Drive/PaoLabShare/FloatingStudy/FAST_SimOut';
%     filename = 'Floating_Hs5Tp15_u10p5_Model.out';
    filename = 'Floating_Hs5Tp8_u15_Model.out';

    file = [dir filesep filename];

    fastout = LoadFastOut(file);
elseif type == 2
    fastout = simout;
else
    error('Make sure your simulation type is defined properly') 
end

fastplot = fastout;

%% Define desired FAST outputs to plot
% outstr = {'Wind1VelX';
%         'BldPitch1';
%              'GenTq';
%             'GenSpeed';
%             'GenPwr';
%             'RotSpeed';
% };

outstr = {'PtfmPitch';
            'TwrBsMyt';
             'LSSTipMys';
             'BldPitch1';
            'RootMyc1';
            'RotSpeed';
            'RotThrust';
    };


%% Plot Outputs
for j = 1:length(outstr)
    fid = figure(400+j);
    f = myplot(fastplot.Time, fastplot.(outstr{j}));
    set(f,'name',outstr{j})
    xlabel('Time (s)')
    title(outstr{j})
    hold on
    
    % Save figures
    runtype = 'VSPI_NoFilt_step_fixed';
    if sfig == 1
        root = '/Volumes/GoogleDrive/My Drive/FOWT_ResearchCollab/Figures/GenericControl_v0/';        
        saveas(fid, [root, outstr{j}, '_', runtype], 'png')
    end
    % Save Data
    if sdata == 1
        root = '/Volumes/GoogleDrive/My Drive/FOWT_ResearchCollab/simout/GenericControl_v0/';
        save([root,runtype], 'fastout')
    end
end



%% Plot using openfast archive script...

% FASTfiles = {'DTU_10MW_RWT.out'};
% FASTfilesDesc = {'Case4'};
% ReferenceFile = 1;
% Channels = {'BldPitch1', 'Azimuth'}
% ShowLegend = 1;
% CustomHdr = [];
% PlotPSDs = 1;
% OnePlot = 0;
% 
% 
% [outData]=PlotFASToutput(FASTfiles,FASTfilesDesc,ReferenceFile,Channels,ShowLegend,CustomHdr,PlotPSDs,OnePlot)
