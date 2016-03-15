function [ TTilde_vec, CV_vec, PTilde_vec, KappaTilde_vec ] ...
    = SimulatedUnitarity( input_args )
%SimulatedUnitarity Simulated quantities for the EOS at unitarity
% The backbone of this program is the theory data for the normalized
% specific heat vs. T/TF.
%

%% include path for virial expansion and ideal EOS
addpath('/Users/Julian/Documents/MIT/MatlabPrograms/VirialExpansion')
addpath('/Users/Julian/Documents/MIT/MatlabPrograms/IdealFermiEOS')

%% Absolute temperature
TAbs = 500*10^(-9); % absolute temperature in Kelvin

%% Bertsch parameter
xi = 0.37;

%% Critical temperature
TcTilde = 0.17;

%% Ideal EOS
[ KappaTilde_ideal, PTilde_ideal, TTilde_ideal, CV_NkB_ideal, beta_mu_vec_ideal ,Z_vec_ideal ] = IdealFermiEOS();

%% load theory data for normalized CV curve
%load('CVTheory.mat');
CV_vec = CV_NkB_ideal;

%% T/TF vector
% [Jump_Value,TC_index] = max(abs(diff(CV_vec)));
% CV_size = size(CV_vec,2);
% TheorySlope = TcTilde/TC_index;
% TTilde_vec = 0:TheorySlope:(CV_size-1)*TheorySlope;

TTilde_vec = TTilde_ideal;

%% normalized Pressure as a function of T/TF
PTilde_vec = 5/3 * cumtrapz(TTilde_vec,CV_vec) + PTilde_ideal(1);

%% normalized Compressibility as a function of T/TF
KappaTilde_vec = 1./(PTilde_vec - 2*TTilde_vec.*CV_vec/3);

%% Beta * Mu as a function of T/TF (don't trust the very low temperature results)
% get initial value for BetaMu in the virial regime

beta_mu = zeros(1,length(TTilde_vec));
beta_mu(1) = beta_mu_vec_ideal(1);
core = (1./(KappaTilde_vec.*(TTilde_vec).^2));
DT = diff(TTilde_vec);
    for i=1:length(TTilde_vec)-1
    beta_mu(length(TTilde_vec)-(i-1)) = -  DT(i:end)*core(i+1:end)' + beta_mu(1);
    end


%% mu/ EF as a function of T/TF (don't trust the very low temperature results)
mu_EF = TTilde_vec.* beta_mu;

%% plot the results

figure(1)
subplot(2,2,1);
plot(TTilde_vec,CV_vec)
xlabel ('$T/T_\mathrm{F}$','interpreter','latex','FontSize',16)
ylabel ('$C_\mathrm{V}/N k_\mathrm{B}$','interpreter','latex','FontSize',16)
xlim([0 1.2])

subplot(2,2,2);
plot(PTilde_vec,KappaTilde_vec)
xlabel ('$P/P_0$','interpreter','latex','FontSize',16)
ylabel ('$\kappa/\kappa_0$','interpreter','latex','FontSize',16)
hold on
plot(PTilde_ideal,KappaTilde_ideal)
hold off
xlim([0 5])

subplot(2,2,3);
plot(TTilde_vec,mu_EF)
xlabel ('$T/T_\mathrm{F}$','interpreter','latex','FontSize',16)
ylabel ('$\mu/E_\mathrm{F}$','interpreter','latex','FontSize',16)
xlim([0.1 6])
hold on
plot(TTilde_ideal,TTilde_ideal.*beta_mu_vec_ideal)
hold off

subplot(2,2,4);
%plot(PTilde_vec,KappaTilde_vec)
xlabel ('$\mu \beta$','interpreter','latex','FontSize',16)
ylabel ('$n(\mu,T)/n_0(\mu,T)$','interpreter','latex','FontSize',16)

figure(2)
plot(TTilde_vec,beta_mu)
xlim([0.1 max(TTilde_vec)])
hold on
plot(TTilde_ideal,beta_mu_vec_ideal)
hold off


figure(3)
plot(TTilde_vec,mu_EF)
xlabel ('$T/T_\mathrm{F}$','interpreter','latex','FontSize',16)
ylabel ('$\mu/E_\mathrm{F}$','interpreter','latex','FontSize',16)
xlim([0 6])
hold on
plot(TTilde_ideal,TTilde_ideal.*beta_mu_vec_ideal)
hold off

figure(4)
plot(TTilde_vec,mu_EF)
xlabel ('$T/T_\mathrm{F}$','interpreter','latex','FontSize',16)
ylabel ('$\mu/E_\mathrm{F}$','interpreter','latex','FontSize',16)
xlim([0.0 0.5])
hold on
plot(TTilde_ideal,TTilde_ideal.*beta_mu_vec_ideal)
hold off

figure(6)
plot(TTilde_ideal,beta_mu_vec_ideal-(beta_mu))
ylim([-0.5 0.5])

end

