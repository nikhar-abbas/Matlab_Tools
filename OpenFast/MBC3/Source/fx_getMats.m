function [matData, data] = fx_getMats(FileNames)
% fx_getMats(FileNames)
% Written by J. Jonkman, NREL
% 19-July-2016: Updated by B. Jonkman (NREL) to convert FAST v8.16 
% linearization files into format expected by mbc3.m
% 22-Jan-2018: Updated by B. Jonkman, (Envision Energy) for BeamDyn linearization
%              converted to a function with data types.
% 
% This m-file is used to read in the data written to multiple FAST linear output
%  (.lin) files, compute the state matrix, [A], at each of the equally-spaced
%  azimuth steps and their azimuth-average, along with their eigenvalues and
%  eigenvectors.
%
% inputs: 
%   FileNames - cell string containing names of FAST linear output files
% outputs:
%   matData   - structure containing computations of FAST linear data.
%   data      - raw data read from the FAST linearization files.
%
% ASSUMPTIONS:
% - all files in FileNames contain the same data structures (i.e., they
%   have the same numbers of inputs, outputs, and states; if a state matrix
%   exists in one, it exists in all)
% - all files in FileNames have the same rotor speed, but have different
%   azimuth angles
% - BeamDyn blade nodes are discretized in the same way for each blade.
% - the states for each module are ordered by module; in each module all of  
%   the displacement states are followed by all the velocity states.
% - descriptions of inputs, outputs, and (BeamDyn) states are triplets if
%   they match in all characters except the blade number. (see
%   findBladeTriplets.m for details)

if nargin < 1 || isempty(FileNames)
    FileNames = {'Test18.1.lin','Test18.2.lin'};
 elseif ~iscell(FileNames)
     FileNames = {FileNames}; % convert (hopefully a) string to cell; it does not...
     d=dir('*.lin');
     for ix=1:length(d)
         FileNames{ix}=d(ix).name;
     end
end 

% Input data from linearization files:
matData.NAzimStep       = length(FileNames);
data(matData.NAzimStep) = ReadFASTLinear(FileNames{matData.NAzimStep}); %we'll read this twice so we can allocate space first; putting it at matData.NAzimStep saves some reallocation later
matData.NumStates       = data(matData.NAzimStep).n_x;
matData.ndof            = data(matData.NAzimStep).n_xx / 2 + data(matData.NAzimStep).n_x - data(matData.NAzimStep).n_xx; %new calc for dof
matData.n_xx            = data(matData.NAzimStep).n_xx; %number of second order states; consider deleting 12/17/18
matData.NumInputs       = data(matData.NAzimStep).n_u;
matData.NumOutputs      = data(matData.NAzimStep).n_y;


%% .................................
% allocate space for these variables
% ..................................
matData.Azimuth   = zeros(matData.NAzimStep, 1);
matData.Omega     = zeros(matData.NAzimStep, 1);
matData.OmegaDot  = zeros(matData.NAzimStep, 1);

if ( matData.NumStates > 0 )
    matData.DescStates = data(matData.NAzimStep).x_desc;
% %     matData.xdop       = zeros(matData.NumEDStates, matData.NAzimStep); consider deleting
% %     matData.xop        = zeros(matData.NumEDStates,                    matData.NAzimStep); consider deleting
% %     matData.A          = zeros(matData.NumEDStates, matData.NumEDStates, matData.NAzimStep); consider deleting  
end

if ( matData.NumInputs > 0 )
    matData.DescCntrlInpt = data(matData.NAzimStep).u_desc;    
    [NumEDInputs, NumHDInputs]  =  getInputDetails(matData);
    matData.NumEDInputs = NumEDInputs;
    matData.NumHDInputs = NumHDInputs;
    if (matData.NumStates>0) 
        matData.B = zeros(matData.NumStates, matData.NumInputs, matData.NAzimStep);
    end
end

if ( matData.NumOutputs > 0 )
    matData.DescOtpt = data(matData.NAzimStep).y_desc;    
    [NumEDOutputs, NumHDOutputs]  =  getOutputDetails(matData);
    matData.NumEDOutputs = NumEDOutputs;
    matData.NumHDOutputs = NumHDOutputs;
end

if ( matData.NumOutputs > 0 )
    matData.DescOutput    = data(matData.NAzimStep).y_desc;    
    if ( matData.NumStates > 0 )
        matData.C         = zeros(matData.NumOutputs, matData.NumStates, matData.NAzimStep);
    end
    if ( matData.NumInputs > 0 )
        matData.D         = zeros(matData.NumOutputs, matData.NumInputs, matData.NAzimStep);
    end
end

%% Reorder state matrices so that all the modules' displacements are first,
%  followed by all the modules' velocities (because mbc assumes that the 
%  first ndof values are velocities).
% % % if ( any(data(matData.NAzimStep).x_rotFrame) > 0 )
% %     [StateOrderingIndx, checkEDstates,]    = getStateOrderingIndx(matData); %Orders xMat_ED Matrix
    [StateOrderingIndx_ED, checkEDstates, checkBDstates, NumEDStates, NumBDStates, NumHDStates, NumExctnStates]    = getStateOrderingIndx(matData); %Save old way
    x_rotFrame         = data(matData.NAzimStep).x_rotFrame(StateOrderingIndx_ED); % Updated by NJ, was not ordering correctly before
    matData.DescStates = data(matData.NAzimStep).x_desc(StateOrderingIndx_ED);% Updated by NJ, was not ordering correctly before 
    if NumHDStates > 0 
        matData.DescStates_HD = data(matData.NAzimStep).x_desc(NumEDStates+1:matData.ndof);
        matData.DescStates = vertcat(matData.DescStates, matData.DescStates_HD);
    end
    matData.NumEDStates = NumEDStates; %added by NJ
    matData.NumEDStates_1 = NumEDStates/2; % Number of first-order ED states
    matData.NumBDStates = NumBDStates; %added by NJ
    matData.NumHDStates = NumHDStates; %added by NJ
    matData.NumExctnStates = NumExctnStates; %added by NJ
    matData.NumRdtnStates = NumHDStates - NumExctnStates; %added by NJ
% % % end

%% .................................
% get data into variables expected by GetMats (concatenate data from
% different azimuths into matrices)
% ..................................

for iFile = 1:matData.NAzimStep

    data(iFile) = ReadFASTLinear(FileNames{iFile}); %   data(iFile) = ReadFASTLinear(FileNames{iFile});
    
    matData.Omega(iFile)   = data(iFile).RotSpeed;
    matData.Azimuth(iFile) = data(iFile).Azimuth*180/pi;
    
    %Splitting matData.A etc into ED and HD components for linear algebra combinations in fx_mbc3
    
        if (isfield(data(iFile), 'A'))
            matData.A_ED(StateOrderingIndx_ED,StateOrderingIndx_ED,iFile) = data(iFile).A(StateOrderingIndx_ED,StateOrderingIndx_ED);
        end
        if (isfield(data(iFile), 'B'))
            matData.B_ED(StateOrderingIndx_ED,:,iFile) = data(iFile).B(StateOrderingIndx_ED,:);
        end
        if (isfield(data(iFile), 'C'))
            matData.C_ED(:,StateOrderingIndx_ED,iFile) = data(iFile).C(:,StateOrderingIndx_ED);
        end
        if (isfield(data(iFile), 'D'))
            matData.D_ED(:,:,iFile) = data(iFile).D(:,:);
        end
        if (isfield(data(iFile), 'x_op'))        
            matData.xop(StateOrderingIndx_ED,iFile) = cell2mat(data(iFile).x_op(StateOrderingIndx_ED));
        end
        if (isfield(data(iFile), 'xdot_op'))
            matData.xdop(StateOrderingIndx_ED,iFile) = cell2mat(data(iFile).xdot_op(StateOrderingIndx_ED));
        end
        
    if matData.NumHDStates > 0 %Splitting matData.A etc into ED and HD components for linear algebra combinations in fx_mbc3

        if (isfield(data(iFile), 'A'))
            matData.A_HD(1 : matData.NumHDStates, 1 : matData.NumHDStates, iFile) = data(iFile).A(matData.NumEDStates+1 : matData.NumEDStates+matData.NumHDStates, matData.NumEDStates+1 : matData.NumEDStates + matData.NumHDStates);
        end
        if (isfield(data(iFile), 'A'))
            matData.A_AS(1 : matData.NumEDStates_1, 1 : matData.NumHDStates, iFile) = data(iFile).A(matData.NumEDStates_1 + 1 : matData.NumEDStates, matData.NumEDStates + 1 : matData.NumEDStates + matData.NumHDStates); %coupling terms
        end
        if (isfield(data(iFile), 'B'))
            matData.B_HD(:, :, iFile) = data(iFile).B(matData.NumEDStates + 1 : (matData.NumEDStates + matData.NumHDStates), matData.NumEDInputs + 1 : (matData.NumEDInputs + matData.NumHDInputs) );
        end
        if (isfield(data(iFile), 'C'))
            matData.C_HD(:, 1 : matData.NumHDStates, iFile) = data(iFile).C(:, matData.NumEDStates + 1 : matData.NumEDStates + matData.NumHDStates);
        end
% % %         if (isfield(data(iFile), 'D'))
% % %             matData.D_HD(:, :, iFile) = data(iFile).D(matData.NumEDInputs + 1 : matData.NumEDInputs + matData.NumHDInputs, :);
% % %         end
        if (isfield(data(iFile), 'x_op'))        
            matData.xop_HD(1 : matData.NumHDStates, iFile) = cell2mat(data(iFile).x_op(matData.NumEDStates + 1 : matData.NumEDStates + matData.NumHDStates));
            matData.xop=vertcat(matData.xop, matData.xop_HD);
        end
        if (isfield(data(iFile), 'xdot_op'))
            matData.xdop_HD(1 : matData.NumHDStates, iFile) = cell2mat(data(iFile).xdot_op(matData.NumEDStates + 1 : matData.NumEDStates + matData.NumHDStates));
            matData.xdop=vertcat(matData.xdop,matData.xop_HD);
        end
    end    
end


for i=1:matData.NumEDStates_1
    col = strfind(matData.DescStates{i},'DOF_GeAz'); % find the starting index of the string 'DOF_GeAz'
    if ( ~isempty(col) )     % true if the matData.DescStates contains the string 'DOF_GeAz'
        matData.Omega(:)    = matData.xdop(i,:)';
        matData.OmegaDot(:) = matData.xdop(i+matData.NumEDStates_1,:)';
        break;
    end
end
for i=1:matData.NumEDStates_1
    col = strfind(matData.DescStates{i},'DOF_DrTr'); % find the starting index of the string 'DOF_DrTr'
    if ( ~isempty(col) )     % true if the matData.DescStates contains the string 'DOF_GeAz'
        matData.Omega(:)    = matData.Omega(:)    + matData.xdop(i,:)'; %This always comes after DOF_GeAz so let's just add it here (it won't get written over later).
        matData.OmegaDot(:) = matData.OmegaDot(:) + matData.xdop(i+matData.NumEDStates_1,:)';
        break;
    end
end

% ----------- Find multi-blade coordinate (MBC) transformation indices ----

%% Find the indices for, state triplets in the rotating frame
%  (note that we avoid the "first time derivative" states)
if any(x_rotFrame) > 0 && checkEDstates > 0  %%%(matData.NumEDStates_1 > 0)   
        [matData.RotTripletIndicesStates] = findBladeTriplets_EDstate(x_rotFrame(1:matData.NumEDStates_1),matData.DescStates(1:matData.NumEDStates_1) ); % looks like it includes BD states right now... check this
% %     elseif (checkBDstates)
% %         [matData.RotTripletIndicesStates] = findBladeTriplets_EDstate(x_rotFrame(1:matData.ndof),matData.DescStates(1:matData.ndof) ); % looks like it includes BD states right now... check this
else 
        [matData.RotTripletIndicesStates] = findBladeTriplets(        x_rotFrame(1:matData.NumEDStates_1),matData.DescStates(1:matData.NumEDStates_1) );
end

%% Find the indices for control input triplets in the rotating frame:
if (matData.NumInputs > 0)
    [matData.RotTripletIndicesCntrlInpt] = findBladeTriplets(data(1).u_rotFrame,matData.DescCntrlInpt );
end

%% Find the indices for output measurement triplets in the rotating frame:
if (matData.NumOutputs > 0 )
    [matData.RotTripletIndicesOutput] = findBladeTriplets(data(1).y_rotFrame, matData.DescOutput );
end
    
% % %% Find the indices for first order only triplets in the rotating frame:% Added by NJ; use for rotating states to make more general, not HD
% % if ( (matData.NumStates - matData.n_xx) > 0 )
% %     [matData.RotTripletIndicesFstOrdStates] = findBladeTriplets(x_rotFrame(matData.NumEDStates_1+NumBDStates/2+1:matData.NumEDStates_1+NumBDStates/2+NumHDStates), matData.DescStates (matData.NumEDStates_1+NumBDStates/2+1:matData.NumEDStates_1+NumBDStates/2+NumHDStates));
% % end

return;
end

%% Reorder state matrices so that all the module's displacements are first,
%  followed by all the modules' velocities (because mbc assumes that the 
%  first ndof values are velocities).
function [StateOrderingIndx_ED, checkEDstates, checkBDstates, NumEDStates, NumBDStates, NumHDStates, NumExctnStates] = getStateOrderingIndx(matData) %old way

NumState = 0;
NumBDStates = 0; %remove
NumHDStates = 0;
NumExctnStates = 0;
NumEDStates = 0;
checkEDstates = false;
checkBDstates = false; %remove
lastModName = '';
mod_nDOFs   = 0;    % number of DOFs in each module
sum_nDOFs   = 0;    % running total of DOFs
indx_start  = 1;    % starting index of the modules

    for i=1:1:matData.NumStates 
        modName = strtok(matData.DescStates{i}); % name of the module whose states we are looking at

        if strcmpi(modName,'ED')
            NumState = NumState + 1; %REMOVE???
            NumEDStates = NumEDStates + 1;
            checkEDstates = true; %moved to this loop
        end

        if strcmpi(modName,'BD_1')||strcmpi(modName,'BD_2')||strcmpi(modName,'BD_3')
            NumState = NumState + 1;
            NumBDStates = NumBDStates + 1;
            checkBDstates = true;
        end	

        if strcmpi(modName,'HD')
            NumState = NumState + 1;
            NumHDStates = NumHDStates + 1;
            ExctnName = strtok(matData.DescStates{i}, 'HD ');
            if strcmpi(ExctnName(1:5),'Exctn')
                NumExctnStates = NumExctnStates + 1;
            end
        end
    end
    
StateOrderingIndx_ED = (1:NumEDStates)';    

    for i=1:2:NumEDStates % there are an even number of states, so we're going to save some time
        
        modName = strtok(matData.DescStates{i}); % name of the module whose states we are looking at

        if ~strcmp(lastModName,modName)
            % this is the start of a new set of DOFs, so we'll set the
            % "first time derivative" descriptions to an empty string.
            StateOrderingIndx_ED(  indx_start           :(indx_start+mod_nDOFs-1)) = sum_nDOFs +                (1:mod_nDOFs);
            StateOrderingIndx_ED( (indx_start+mod_nDOFs):(i - 1))                  = sum_nDOFs + matData.ndof + (1:mod_nDOFs);

            % reset for a new module
            sum_nDOFs = sum_nDOFs + mod_nDOFs;
            mod_nDOFs = 0;
            indx_start = i;
            lastModName = modName;
        end
        mod_nDOFs = mod_nDOFs+1;
        
    end
    
    StateOrderingIndx_ED(  indx_start           :(indx_start+mod_nDOFs-1)) = sum_nDOFs +                (1:mod_nDOFs);
    StateOrderingIndx_ED( (indx_start+mod_nDOFs):NumEDStates)        = sum_nDOFs + NumEDStates/2 + (1:mod_nDOFs);
    
        % ED has the blade number in the state description twice, so we
        % have to check the strings differently. We'll note that here:
% remove % %     if strcmpi(lastModName,'ED')
% % %         checkEDstates = true;
% % %     else
% % %         checkEDstates = false;
% remove % %     end

end

function [NumEDInputs, NumHDInputs] = getInputDetails(matData) 

NumEDInputs = 0;
NumHDInputs = 0;
p=length(matData.DescCntrlInpt );
    for i=1:1:p
        modName = strtok(matData.DescCntrlInpt{i}); % name of the module whose states we are looking at

        if strcmpi(modName,'IfW')
            NumEDInputs = NumEDInputs + 1;
        end
        
        if strcmpi(modName,'ED')
            NumEDInputs = NumEDInputs + 1;
        end
        
        if strcmpi(modName,'HD')
            NumHDInputs = NumHDInputs + 1;
        end
    end
    
end 

function [NumEDOutputs, NumHDOutputs] = getOutputDetails(matData) 

NumEDOutputs = 0;
NumHDOutputs = 0;
p=length(matData.DescOtpt );
    for i=1:1:p
        modName = strtok(matData.DescOtpt{i}); % name of the module whose states we are looking at

        if strcmpi(modName,'SrvD')
            NumEDOutputs = NumEDOutputs + 1;
        end
        
        if strcmpi(modName,'ED')
            NumEDOutputs = NumEDOutputs + 1;
        end
        
        if strcmpi(modName,'HD')
            NumHDOutputs = NumHDOutputs + 1;
        end
    end
    
end 