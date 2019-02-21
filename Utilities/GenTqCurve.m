% Create Generator Torque Curve Plot

%% NREL 5MW OC3-Hywind
VS_RtGnSp 	  = 	121.6805;
VS_RtPwr      = 	5296610.0;
PC_RefSpd     =     122.9096;
VS_CtInSp     =     70.16224;
VS_Rgn2Sp     =     91.21091;
VS_Rgn2K      =     2.332287;
VS_SlPc       =     10.0;


VS_SySp    = VS_RtGnSp/ (1.0 +  0.01*VS_SlPc);
VS_Slope15 = (VS_Rgn2K*VS_Rgn2Sp*VS_Rgn2Sp)/(VS_Rgn2Sp - VS_CtInSp);
VS_Slope25 = (VS_RtPwr/PC_RefSpd)/(VS_RtGnSp - VS_SySp);
VS_TrGnSp = ( VS_Slope25 - sqrt(VS_Slope25*(VS_Slope25 - 4.0*VS_Rgn2K*VS_SySp)))/(2.0*VS_Rgn2K);


%% DTU 10 MW
VS_RtGnSp 	  = 	42.411;             % rpm
VS_RtPwr      = 	1.0E+07	;           % W
PC_RefSpd     =     50.26548;           % rpm        
VS_CtInSp     =     31.415;             % rpm
VS_Rgn2Sp     =     42.411;             % rpm
VS_Rgn2K      =     104.1;              % N-m/(rpm)^2
VS_SlPc       =     10.0;


VS_SySp    = VS_RtGnSp/ (1.0 +  0.01*VS_SlPc);
VS_Slope15 = (VS_Rgn2K*VS_Rgn2Sp*VS_Rgn2Sp)/(VS_Rgn2Sp - VS_CtInSp);
VS_Slope25 = (VS_RtPwr/PC_RefSpd)/(VS_RtGnSp - VS_SySp);
VS_TrGnSp = ( VS_Slope25 - sqrt(VS_Slope25*(VS_Slope25 - 4.0*VS_Rgn2K*VS_SySp)))/(2.0*VS_Rgn2K);
%%
Om = 0:.1:150;

% for i = 1:length(Om)
%     if Om(i) < Om_ci
%         Tau_g(i) = 0;
%     elseif Om(i) <= Om_r
%         Tau_g(i) = K*Om(i)^2;
%     else
%         Tau_g(i) = Tau_r;
%     end
% end

for i = 1:length(Om)
    if Om(i) >= VS_RtGnSp 
        GenTrq(i) = VS_RtPwr/PC_RefSpd;
    elseif Om(i) <= VS_CtInSp                   % We are in region 1 - torque is zero
         GenTrq(i) = 0.0;
    elseif Om(i) <  VS_Rgn2Sp                   % We are in region 1 1/2 - linear ramp in torque from zero to optimal
         GenTrq(i) = VS_Slope15*( Om(i) - VS_CtInSp );
    elseif Om(i) <  VS_TrGnSp                 % We are in region 2 - optimal torque is proportional to the square of the generator speed
         GenTrq(i) = VS_Rgn2K*Om(i)*Om(i);
    else                                                     % We are in region 2 1/2 - simple induction generator transition region
         GenTrq(i) = VS_Slope25 * (Om(i) - VS_SySp);
    end
end
 
%% Plotting
close all
fig = myplot(Om,GenTrq);
% fig = AIAA_pplot(Om*30/pi,GenTrq./1000); hold on
% fig = AIAA_pplot(Om*30/pi, VS_Rgn2K*Om.^2./1000,'k--');
% Leg = legend('Variable Speed Controller','Optimal','location','northwest');
Leg.FontSize = 14;
ylabel('Generator Torque, Nm')
xlabel('Generator Speed, rpm')
title('Generator Torque Curve','fontsize',16)


%% Save plot
% figdir = '~/Documents/Presentations/AIAA_2019/Figures';
% figname = 'GenTqCv.eps';
% saveas(fig,[figdir filesep figname],'epsc')
