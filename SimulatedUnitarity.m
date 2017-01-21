function [ TTilde_vec, CV_vec, PTilde_vec, KappaTilde_vec, mu_EF_vec, BetaMu_vec ] ...
    = SimulatedUnitarity( );
%SimulatedUnitarity Simulated quantities for the EOS at unitarity
% The backbone of this program is the theory data for the normalized
% specific heat vs. T/TF.
%

%% include path for virial expansion and ideal EOS
addpath('/Users/RanchoP/Dropbox (MIT)/Github/VirialExpansion')
addpath('/Users/RanchoP/Dropbox (MIT)/Github/IdealFermiEOS')

%% Bertsch parameter
xi = 0.37;

%% Critical temperature
TcTilde = 0.17;

%% load theory data for normalized CV curve
load('CVTheory_extended.mat');
%CV_vec = CV_vec';
CV_vec = CV_vec(1:10^5);

%% T/TF vector
[Jump_Value,TC_index] = max(abs(diff(CV_vec)));
CV_size = size(CV_vec,2);
TheorySlope = TcTilde/TC_index;
TTilde_vec = 0:TheorySlope:(CV_size-1)*TheorySlope;

%% normalized Pressure as a function of T/TF
PTilde_vec = 5/3 * cumtrapz(TTilde_vec,CV_vec) + xi;

%% normalized Compressibility as a function of T/TF
KappaTilde_vec = 1./(PTilde_vec - 2*TTilde_vec.*CV_vec/3);

%% Ideal EOS
[ KappaTilde_ideal, PTilde_ideal, TTilde_ideal, CV_NkB_ideal,  beta_mu_vec_ideal ,Z_vec_ideal ] = IdealFermiEOS( );

%% Beta * Mu as a function of T/TF (don't trust the very low temperature results)
% get initial value for BetaMu in the virial regime
[ KappaTilde_Virial, PTilde_Virial, TTilde_Virial, CI_NkF_Virial, BetaMu_vec_Virial, Z_vec_Virial ] = ...
    VirialUnitarity('LogPoints',10^7);

[Diviation,Virial_index] = min(abs(max(TTilde_vec)-TTilde_Virial));

BetaMu_initial = BetaMu_vec_Virial(Virial_index);

BetaMu_vec = zeros(1,length(TTilde_vec));
BetaMu_vec(1) = BetaMu_initial;
Int_core = (1./(KappaTilde_vec.*(TTilde_vec).^2));
DT = diff(TTilde_vec);
    for i=1:length(TTilde_vec)-1
    BetaMu_vec(i+1) = +  DT(i:end)*Int_core(i+1:end)' + BetaMu_vec(1);
    end


%% mu/ EF as a function of T/TF (don't trust the very low temperature results)
mu_EF_vec = TTilde_vec.* BetaMu_vec;

%% plot the results

figure(1)
subplot(3,1,1);
plot(TTilde_vec,CV_vec,'k')
xlabel ('$T/T_\mathrm{F}$','interpreter','latex','FontSize',16)
ylabel ('$C_\mathrm{V}/N k_\mathrm{B}$','interpreter','latex','FontSize',16)
hold on
plot(TTilde_ideal,CV_NkB_ideal,'r')
hold off
xlim([0 1.2])

subplot(3,1,2);
plot(PTilde_vec,KappaTilde_vec,'k')
xlabel ('$P/P_0$','interpreter','latex','FontSize',16)
ylabel ('$\kappa/\kappa_0$','interpreter','latex','FontSize',16)
hold on
plot(PTilde_Virial,KappaTilde_Virial,'g')
plot(PTilde_ideal,KappaTilde_ideal,'r')
hold off
xlim([0 5])

subplot(3,1,3);
plot(TTilde_vec,mu_EF_vec,'k')
xlabel ('$T/T_\mathrm{F}$','interpreter','latex','FontSize',16)
ylabel ('$\mu/E_\mathrm{F}$','interpreter','latex','FontSize',16)
xlim([0.1 1.6])
hold on
plot(TTilde_Virial,TTilde_Virial.*BetaMu_vec_Virial,'g')
plot(TTilde_ideal,TTilde_ideal.*beta_mu_vec_ideal,'r')
hold off

%subplot(2,2,4);
%plot(PTilde_vec,KappaTilde_vec)
%xlabel ('$\mu \beta$','interpreter','latex','FontSize',16)
%ylabel ('$n(\mu,T)/n_0(\mu,T)$','interpreter','latex','FontSize',16)

% figure(2)
% plot(TTilde_vec,BetaMu_vec)
% xlim([0.1 max(TTilde_vec)])
% hold on
% plot(TTilde_Virial,BetaMu_vec_Virial)
% plot(TTilde_ideal,beta_mu_vec_ideal)
% hold off
% 
% figure(3)
% plot(TTilde_vec,mu_EF)
% xlabel ('$T/T_\mathrm{F}$','interpreter','latex','FontSize',16)
% ylabel ('$\mu/E_\mathrm{F}$','interpreter','latex','FontSize',16)
% xlim([0.0 6])
% hold on
% plot(TTilde_ideal,TTilde_ideal.*beta_mu_vec_ideal)
% hold off



end

