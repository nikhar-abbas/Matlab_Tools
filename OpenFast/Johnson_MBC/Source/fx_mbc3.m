function [MBC, matData, FAST_linData] = fx_mbc3( FileNames ) 
% MBC: Multi-Blade Coordinate Transformation for a turbine with 3-blade rotor
%
% Developed by Gunjit Bir, NREL (303-384-6953, gunjit_bir@nrel.gov)
%
% Last update: 08/30/2006
%    29-Jan-2018  - B. Jonkman, Envision Energy & J. Jonkman, NREL
%----------------------------------------------------------------------------------------
%
% Objectives:
% 1. Given state-space matrices (A,B) and output matrices (C,D), defined partly in the
%    rotoating frame and partly in the fixed frame, transform these matrices to the fixed
%    coordinate frame using multi-blade coordinate transformation (MBC). The transformned
%    matrices are MBC.A, MBC.B, MBC.C, and MBC.D.
%
% 2. Given second-order matrices (M,Dmp,K), control input matrix (F), 
%    displacement output matrix (Dc), and velocity output matrix (Vc), transform
%    these to the fixed coordinate frame using multi-blade coordinate transformation (MBC).
%    The transformed matrices are MBC.M, MBC.Dmp, MBC.K, MBC.F, MBC.Dc, and MBC.Vc.
%
% 3. Azimuth-average the MBC.A matrix and compute its eigenvalues and eigenvectors.  The
%    eigenvectors, referred to the fixed coordinates, are presented in both complex and
%    amplitude-phase forms. The eigenvalues, also referred to the fixed coordinates, are
%    presented in complex and damping-ratio/damped-frequency forms.
%
% ***Disclaimer: This code is still in the developmental stage and no guarantee is given
%    as to its proper functioning or accuracy of its results.
%
% ------------ INPUTS   (these are fields in matData [see fx_getMatx.m]) ---
% NumInputs  : total num of control inputs
% NumOutputs : total num of outputs
% NumStates  : number of states
%
% RotTripletIndicesStates   : State triplets in rotating frame (matrix of size rotating_dof_types*3)
% RotTripletIndicesCntrlInpt: Control-input triplets in rotating frame (matrix of size rotating_control_input_types*3)
% RotTripletIndicesOutput   : Output triplets in rotating frame (matrix of size rotating_output_types*3)
%
% DescStates : description of states associated with input matrices (FAST-specific) %%
%
% A, B, C, D:                     1st-order input matrices
% MassMat, DampMat, StffMat, FMat, VelCMat, DspCMat: 2nd-order input matrices
% Omega     : Vector of rotor speeds at specified azimuths (rad/sec)
% OmegaDot  : Vector of rotor accelerations at specified azimuths (rad/sec2)
% ndof      : total number of degrees of freedom
% NumStates : number of states
% Azimuth   : vector of azimuth positions (in deg)
%
% --------------------- OUTPUTS ----------------------------------------------------
% MBC.A,MBC.B         : state-space matrices transformed to the fixed frame
% MBC.C,MBC.D         : output matrices transformed to the fixed frame
% MBC.M,MBC.Dmp,MBC.K : second-order mass, damping/gyroscopic and stiffness matrices transformed to the fixed frame
% MBC.F               : control input matrix transformed to the fixed frame
% MBC.Dc,MBC.Vc       : displacement and velocity output matrices (Vc) transformed to the fixed frame
% -----------------------------------------------------------------------------------------
%**NOTE 1: All inputs and output matrices may not be present.  For example, user may supply or may be interested
%          in multi-blade coordinate transformation of A matrix only.  In this case, only MBC.A matrix along with
%          fixed-frame eigensolution will be genertaed as outputs.  The code checks for consistency and completness
%          of selected inputs and generates associated outputs.
%
% Equation numbers are from this document: https://nwtc.nrel.gov/system/files/MBC3.pdf
% -----------------------------------------------------------------------------------------

fprintf( '\n  Running %s\n\n', 'mbc3 (v2.0, 29-Jan-2018)' );

[matData, FAST_linData] = fx_getMats(FileNames);

MBC.DescStates = matData.DescStates; % save this in the MBC type for possible campbell_diagram processing later 

%% ---------- Multi-Blade-Coordinate transformation -------------------------------------------
if ~isfield(matData,'RotTripletIndicesStates')
    error('*** There are no rotating states. MBC transformation, therefore, cannot be performed.');
end

[n_RotTripletStates,nb] = size(matData.RotTripletIndicesStates);
if(nb ~= 3)
    error('**ERROR: the number of column vectors in matData.RotTripletIndicesStates must equal 3, the num of blades');
elseif(n_RotTripletStates*nb > matData.ndof)
    error('**ERROR: the rotating dof exceeds the total num of dof');
end

     new_seq_dof_ED    = get_new_seq(matData.RotTripletIndicesStates,matData.NumEDStates_1); % these are the first ndof states (not "first time derivative" states); these values are used to calculate matrix transformations
     new_seq_states_ED = [new_seq_dof_ED  new_seq_dof_ED+matData.NumEDStates_1]; % add the remaining ones (assumes ordering of displacements and velocities in state matrices); these values are used to calculate matrix transformations 

     if matData.NumHDStates > 1
         matData.NumStates = matData.ndof + matData.NumEDStates_1;
%          new_seq_dof_Exctn    = get_new_seq(0,matData.NumExctnStates/2);
%          % these are the first ndof states (not "first time derivative"
%          states); <------------------------ NOT CURRENTLY USED TO
%          CALCULATE MATRICIES
%          new_seq_states_Exctn = [new_seq_dof_Exctn  new_seq_dof_Exctn+matData.NumEDStates_1]; % add the remaining ones (assumes ordering of displacements and velocities in state matrices)
%          new_seq_dof_Rdtn    = get_new_seq(0,matData.NumEDStates_1); % these are the first ndof states (not "first time derivative" states)
%          new_seq_states_Rdtn = [new_seq_dof_Rdtn  new_seq_dof_Rdtn+matData.NumEDStates_1]; % add the remaining ones (assumes ordering of displacements and velocities in state matrices)
     else
         matData.NumStates = matData.NumEDStates;
     end
     
% end

if isfield(matData,'RotTripletIndicesCntrlInpt')
    [n_RotTripletInputs,nb] = size(matData.RotTripletIndicesCntrlInpt);
    if(nb ~= 3)
        error('**ERROR: the number of column vectors in RotTripletIndicesCntrlInpt must equal 3, the num of blades');
    end
    new_seq_inp_ED = get_new_seq(matData.RotTripletIndicesCntrlInpt,matData.NumEDInputs);
else
    n_RotTripletInputs = 0; % number of rotating-frame control triplets
    new_seq_inp_ED = 1:matData.NumEDInputs;
end

if isfield(matData,'RotTripletIndicesOutput')
    [n_RotTripletOutputs,nb] = size(matData.RotTripletIndicesOutput);
    if(nb ~= 3)
        error('**ERROR: the number of column vectors in RotTripletIndicesOutput must equal 3, the num of blades');
    end
    new_seq_out = get_new_seq(matData.RotTripletIndicesOutput,matData.NumOutputs);
else
    n_RotTripletOutputs = 0; % number of rotating-frame output triplets
    new_seq_out = 1:matData.NumOutputs;
end


if isfield(matData,'RotTripletIndicesHydroStates') %added by NJ no rotating HD states
    [n_RotTripletHydroStates,nb] = size(matData.RotTripletIndicesHydroStates);
    if(nb ~= 3)
        error('**ERROR: the number of column vectors in RotTripletIndicesHydroStates must equal 3, the num of blades');
    end
    new_seq_HD_dof = get_new_seq(matData.RotTripletIndicesHydroStates,matData.NumHydroStates);
else
    n_RotTripletHydroStates = 0; % number of rotating-frame output triplets
    new_seq_ = 1:matData.NumHDStates;
end



n_FixFrameStates_ED  = matData.NumEDStates_1       - n_RotTripletStates*nb;  % fixed-frame dof updated by NJ
n_FixFrameInputs_ED  = matData.NumInputs  - n_RotTripletInputs*nb  -  matData.NumHDInputs;  % fixed-frame control inputs

n_FixFrameOutputs = matData.NumOutputs - n_RotTripletOutputs*nb; % fixed-frame outputs

if isfield(matData,'A_ED')
    MBC.AvgA = zeros(matData.NumStates); % initalize matrix
end

if ( size(matData.Omega) ~= matData.NAzimStep)
   error('**ERROR: the size of Omega vector must equal matData.NAzimStep, the num of azimuth steps');
end
if ( size(matData.OmegaDot) ~= matData.NAzimStep)
   error('**ERROR: the size of OmegaDot vector must equal matData.NAzimStep, the num of azimuth steps');
end

% begin azimuth loop 
for iaz = matData.NAzimStep:-1:1  
    T1q = 0;
    T2q = 0;
    T1qv = 0;

    %(loop backwards so we don't reallocate memory each time [i.e. variables with iaz index aren't getting larger each time])

    % compute azimuth positions of blades:
    az = matData.Azimuth(iaz)*pi/180.0 + 2*pi/nb* (0:(nb-1)) ; % Eq. 1, azimuth in radians

    % get rotor speed squared
    OmegaSquared = matData.Omega(iaz)^2;

    % compute transformation matrices
    cos_col = cos(az(:));
    sin_col = sin(az(:));
    
    tt = [ones(3,1), cos_col, sin_col];         % Eq. 9, t_tilde
    ttv = get_tt_inverse(sin_col, cos_col);     % inverse of tt (computed analytically in function below)
    
    %---
    T1 = eye(n_FixFrameStates_ED);                 % Eq. 11
    for ii = 1:n_RotTripletStates
        T1 = blkdiag(T1, tt);
    end

    T1v = eye(n_FixFrameStates_ED);                % inverse of T1
    for ii = 1:n_RotTripletStates
        T1v = blkdiag(T1v, ttv);
    end

    T2 = zeros(n_FixFrameStates_ED);               % Eq. 14
    tt2 = [zeros(3,1), -sin_col,  cos_col];     % Eq. 16 a
    for ii = 1:n_RotTripletStates
        T2 = blkdiag(T2, tt2);
    end

    T3 = zeros(n_FixFrameStates_ED);               % Eq. 15
    tt3 = [zeros(3,1), -cos_col, -sin_col];     % Eq. 16 b
    for ii = 1:n_RotTripletStates
        T3 = blkdiag(T3, tt3);
    end
    
    %---
    T1c = eye(n_FixFrameInputs_ED);                % Eq. 21
    for ii = 1:n_RotTripletInputs;
        T1c = blkdiag(T1c, tt);
    end

    T1ov = eye(n_FixFrameOutputs);              % inverse of Tlo (Eq. 23)
    for ii = 1:n_RotTripletOutputs
        T1ov = blkdiag(T1ov, ttv);
    end

    if matData.NumHDStates>0
        T1q = eye(matData.NumHDStates);                %Added by NJ 9/13/18;
        T2q = zeros(matData.NumHDStates);                %Added by NJ 9/13/18;
        T1qv = eye(matData.NumHDStates);                % inverse of T1q
        T1qc = eye(matData.NumHDInputs);                % inverse of T1q
    end
    
MBC.A = zeros(matData.NumStates); % initalize matrix
% mbc transformation of first-order matrices
%  if ( MBC.EqnsOrder == 1 ) % activate later

    if isfield(matData,'A_ED')
            % Eq. 29, assuming
            % xAMat( 1:matData.ndof, 1:matData.NumStates ) = 0 and
            % xAMat( 1:matData.ndof, (matData.ndof+1):matData.NumStates) = I
        xAMat_ED(:,:) = matData.A_ED(new_seq_states_ED,new_seq_states_ED,iaz); %--
        AK = xAMat_ED( (matData.NumEDStates_1 + 1) : matData.NumEDStates,             1 : matData.NumEDStates_1);
        AC = xAMat_ED( (matData.NumEDStates_1 + 1) : matData.NumEDStates, matData.NumEDStates_1 + 1 : matData.NumEDStates);
               
        MBC.A(new_seq_states_ED,new_seq_states_ED,iaz) = ...
               [zeros(matData.NumEDStates_1),   eye(matData.NumEDStates_1);
                T1v*(AK*T1 +   matData.Omega(iaz)*AC*T2 - OmegaSquared*T3 - matData.OmegaDot(iaz)*T2), ...
                T1v*(AC*T1 - 2*matData.Omega(iaz)*T2)];
            
        if matData.NumHDStates > 0 %matData.NumExctnStates && matData.NumExctnStates > 0
            xMat_HD(:,:) =  matData.A_HD(:,:,iaz);
            AQ = xMat_HD; 
            AS = matData.A_AS(:,:,iaz);
            MBC.A(matData.NumEDStates + 1:matData.NumStates,matData.NumEDStates + 1:matData.NumStates, iaz) = ...
			[T1qv*(AQ*T1q - matData.Omega(iaz)*T2q)]; % Trivial for HD only, no rotation
            MBC.A(matData.NumEDStates_1 + 1:matData.NumEDStates, matData.NumEDStates+1:matData.NumStates, iaz) = ...
            [T1v*(AS*T1q)]; % Trival for HD only, no rotation
        end
    end
          

    if isfield(matData,'B_ED')
            % Eq. 30
        xBMat_ED = matData.B_ED(new_seq_states_ED,new_seq_inp_ED,iaz); %--

        B1 = xBMat_ED(1:matData.NumEDStates_1, 1 : matData.NumEDInputs);
        B2 = xBMat_ED(matData.NumEDStates_1+1:matData.NumEDStates, 1 : matData.NumEDInputs);
        MBC.B(new_seq_states_ED,new_seq_inp_ED,iaz) = [T1v*B1; T1v*B2] * T1c;
        
        if matData.NumHDStates > 0
           xBMat_HD = matData.B_HD(:,:,iaz);
            B3 = xBMat_HD; %Added by NJ 8/16/18; xAmat for Hydro states
           MBC.B(matData.NumEDStates + 1:matData.NumStates, matData.NumEDInputs + 1  : (matData.NumEDInputs + matData.NumHDInputs) ,iaz) = [T1qv*B3] * T1qc; % updated for hydro states  
        end
    end

    if isfield(matData,'C_ED')
            % Eq. 31

        MBC.C(new_seq_out, new_seq_states_ED, iaz) = ...
                     T1ov * matData.C_ED(new_seq_out,new_seq_states_ED,iaz) * ...
                     [T1, zeros(matData.NumEDStates_1); matData.Omega(iaz)*T2, T1]; 
        if matData.NumHDStates > 0
          MBC.C(:, matData.NumEDStates + 1:matData.NumStates, iaz) = ...
              T1ov * matData.C_HD(:,:,iaz) * T1qc; % updated for hydro states  
        end
    end

    if isfield(matData,'D_ED')
           % Eq. 32
        MBC.D(new_seq_out,new_seq_inp_ED,iaz) = T1ov * matData.D_ED(new_seq_out,new_seq_inp_ED,iaz) * T1c;

% % %          
% % %             MBC.D((matData.NumEDOutputs + 1 : matData.NumEDOutputs + matData.NumHDOutputs), : ,iaz) = ...
% % %               matData.D_HD((matData.NumEDOutputs + 1:matData.NumHDOutputs),:,iaz); % updated for hydro states  
            MBC.D(:, 7 ,iaz) = ...
              matData.D_ED(:, 7 ,iaz); % NEED TO FIX!!!!!!!!!!!!!!!!!!!!!!              
% % %         end
    end

% mbc transformation of second-order matrices <--- comment out for new mbc3

end   % end of azimuth loop



%------------- Eigensolution and Azimuth Averages -------------------------
if isfield(MBC,'A')
    MBC.AvgA = mean(MBC.A,3); % azimuth-average of azimuth-dependent MBC.A matrices
     MBC.eigSol = eiganalysis(MBC.AvgA,matData.NumEDStates);
end

if isfield(MBC,'B')
    MBC.AvgB = mean(MBC.B,3); % azimuth-average of azimuth-dependent MBC.B matrices
end

if isfield(MBC,'C')
    MBC.AvgC = mean(MBC.C,3); % azimuth-average of azimuth-dependent MBC.C matrices
end

if isfield(MBC,'D')
    MBC.AvgD = mean(MBC.D,3); % azimuth-average of azimuth-dependent MBC.D matrices
end

MBC.NumEDStates = matData.NumEDStates;
% ----------- Clear unneeded variables -------------------------------
  disp('  ');
  disp(' Multi-Blade Coordinate transformation completed ');
%-----------------------------------------------------------
return;
end

%% ------------------------------------------------------------------------
% compute the inverse of tt = [ones(3,1), cos_col, sin_col]
function [ttv] = get_tt_inverse(sin_col, cos_col)

    c1 = cos_col(1);
    c2 = cos_col(2);
    c3 = cos_col(3);
    
    s1 = sin_col(1);
    s2 = sin_col(2);
    s3 = sin_col(3);

    
    ttv = [ c2*s3 - s2*c3,  c3*s1 - s3*c1, c1*s2 - s1*c2
               s2 - s3 ,       s3 - s1,       s1 - s2
               c3 - c2 ,       c1 - c3,       c2 - c1 ] / (1.5*sqrt(3));

    return
end

%% ------------------------------------------------------------------------
% create a sequence where the non-rotating values are first, and are then 
% followed by the rotating series with b1, b2, b3 triplets:
function [new_seq] = get_new_seq(rot_triplet,ntot,NumED)
%  rot_q_triplet is size n x 3
% if rot_triplet == 0
%  d=4;
% else
      non_rotating = true(ntot,1);
      non_rotating(rot_triplet(:)) = false; % if they are rotating, set them false;
      new_seq = [find(non_rotating); reshape( rot_triplet', numel(rot_triplet), 1)];
      return    
%end
    
    
end

