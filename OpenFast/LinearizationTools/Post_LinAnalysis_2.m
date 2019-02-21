%% Post_LinAnalysis_2



% name states
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

    
    
    st{j}(sta-1:end) = [];
    st{j}(1:3) = [];
    st{j} = strcat(st{j},' - ',num2str(j));
end
% for j = 1:length(linout)
    

%% some plots
st2 = flipud(st);
cats = categorical(st,st2);
ctrmag = abs(linavg.SVD.V(:,1));
fig = figure(2)
b = barh(cats,[ctrmag]')
fig.Children.FontSize = 12;
% leg = legend('\tau','\beta','\tau\beta')
% leg.FontSize = 20;

%%
FullState_desc = {'ED Platform horizontal surge translation DOF (internal DOF index = DOF_Sg), m',...
        'ED Platform horizontal sway translation DOF (internal DOF index = DOF_Sw), m',...
        'ED Platform vertical heave translation DOF (internal DOF index = DOF_Hv), m',...
        'ED Platform roll tilt rotation DOF (internal DOF index = DOF_R), rad',...
        'ED Platform pitch tilt rotation DOF (internal DOF index = DOF_P), rad',...
        'ED Platform yaw rotation DOF (internal DOF index = DOF_Y), rad',...
        'ED 1st tower fore-aft bending mode DOF (internal DOF index = DOF_TFA1), m',...
        'ED 1st tower side-to-side bending mode DOF (internal DOF index = DOF_TSS1), m',...
        'ED 2nd tower fore-aft bending mode DOF (internal DOF index = DOF_TFA2), m',...
        'ED 2nd tower side-to-side bending mode DOF (internal DOF index = DOF_TSS2), m',...
        'ED Nacelle yaw DOF (internal DOF index = DOF_Yaw), rad',...
        'ED Variable speed generator DOF (internal DOF index = DOF_GeAz), rad',...
        'ED Drivetrain rotational-flexibility DOF (internal DOF index = DOF_DrTr), rad',...
        'ED 1st flapwise bending-mode DOF of blade 1 (internal DOF index = DOF_BF(1,1)), m',...
        'ED 1st flapwise bending-mode DOF of blade 2 (internal DOF index = DOF_BF(2,1)), m',...
        'ED 1st flapwise bending-mode DOF of blade 3 (internal DOF index = DOF_BF(3,1)), m',...
        'ED 1st edgewise bending-mode DOF of blade 1 (internal DOF index = DOF_BE(1,1)), m',...
        'ED 1st edgewise bending-mode DOF of blade 2 (internal DOF index = DOF_BE(2,1)), m',...
        'ED 1st edgewise bending-mode DOF of blade 3 (internal DOF index = DOF_BE(3,1)), m',...
        'ED 2nd flapwise bending-mode DOF of blade 1 (internal DOF index = DOF_BF(1,2)), m',...
        'ED 2nd flapwise bending-mode DOF of blade 2 (internal DOF index = DOF_BF(2,2)), m',...
        'ED 2nd flapwise bending-mode DOF of blade 3 (internal DOF index = DOF_BF(3,2)), m',...
        'ED First time derivative of Platform horizontal surge translation DOF (internal DOF index = DOF_Sg), m/s',...
        'ED First time derivative of Platform horizontal sway translation DOF (internal DOF index = DOF_Sw), m/s',...
        'ED First time derivative of Platform vertical heave translation DOF (internal DOF index = DOF_Hv), m/s',...
        'ED First time derivative of Platform roll tilt rotation DOF (internal DOF index = DOF_R), rad/s',...
        'ED First time derivative of Platform pitch tilt rotation DOF (internal DOF index = DOF_P), rad/s',...
        'ED First time derivative of Platform yaw rotation DOF (internal DOF index = DOF_Y), rad/s',...
        'ED First time derivative of 1st tower fore-aft bending mode DOF (internal DOF index = DOF_TFA1), m/s',...
        'ED First time derivative of 1st tower side-to-side bending mode DOF (internal DOF index = DOF_TSS1), m/s',...
        'ED First time derivative of 2nd tower fore-aft bending mode DOF (internal DOF index = DOF_TFA2), m/s',...
        'ED First time derivative of 2nd tower side-to-side bending mode DOF (internal DOF index = DOF_TSS2), m/s',...
        'ED First time derivative of Nacelle yaw DOF (internal DOF index = DOF_Yaw), rad/s',...
        'ED First time derivative of Variable speed generator DOF (internal DOF index = DOF_GeAz), rad/s',...
        'ED First time derivative of Drivetrain rotational-flexibility DOF (internal DOF index = DOF_DrTr), rad/s',...
        'ED First time derivative of 1st flapwise bending-mode DOF of blade 1 (internal DOF index = DOF_BF(1,1)), m/s',...
        'ED First time derivative of 1st flapwise bending-mode DOF of blade 2 (internal DOF index = DOF_BF(2,1)), m/s',...
        'ED First time derivative of 1st flapwise bending-mode DOF of blade 3 (internal DOF index = DOF_BF(3,1)), m/s',...
        'ED First time derivative of 1st edgewise bending-mode DOF of blade 1 (internal DOF index = DOF_BE(1,1)), m/s',...
        'ED First time derivative of 1st edgewise bending-mode DOF of blade 2 (internal DOF index = DOF_BE(2,1)), m/s',...
        'ED First time derivative of 1st edgewise bending-mode DOF of blade 3 (internal DOF index = DOF_BE(3,1)), m/s',...
        'ED First time derivative of 2nd flapwise bending-mode DOF of blade 1 (internal DOF index = DOF_BF(1,2)), m/s',...
        'ED First time derivative of 2nd flapwise bending-mode DOF of blade 2 (internal DOF index = DOF_BF(2,2)), m/s',...
        'ED First time derivative of 2nd flapwise bending-mode DOF of blade 3 (internal DOF index = DOF_BF(3,2)), m/s',...
        }

