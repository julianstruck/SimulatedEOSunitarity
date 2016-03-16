[ TTilde_vec, CV_vec, PTilde_vec, KappaTilde_vec, mu_EF_vec ] ...
    = SimulatedUnitarity( );

%% Input parameters

% Select T/TF and density in center of trap
TTilde0 = 0.1;
n0 = 10^(16); % in 1/m^3

% Select axial trapping frequency
omega_y = 2*pi * 23.9;

% Top pixel size on the atoms
px_top = 1.44*10^(-6); % not used currently)

% physical constants
hbar = 1.0545718*10^(-34);
kB = 1.38064852*10^(-23);

% Lithium parameters
mLi = 9.9883414*10^(-27);

% peak fermi energy
EF0 = hbar^2/(2*mLi)*(3*pi^2*n0)^(2/3);

% peak fermi temperature
TF0 = EF0/kB;

% Temperature of the gas
Tabs = TTilde0 * TF0;

%% Find chemical potential for given T/TF and density in center of trap
[diff,SelectIndex] = min(abs(TTilde_vec - TTilde0));
mu_EF_select = mu_EF_vec(SelectIndex);
mu0 = mu_EF_select*EF0;

%% chemical potential as a function of T/TF
mu_TTilde_vec = mu_EF_vec.*(TTilde_vec).^(-1)*kB*Tabs;

%% chemical potential in LDA
select_mu_TTilde_vec = mu_TTilde_vec(mu0 - mu_TTilde_vec > 0);
z_TTilde_vec = sqrt(2*(mu0 - select_mu_TTilde_vec)/(mLi*omega_y^2));

%% density as a function of T/TF
n_vec_TTilde = 1/(3*pi^2) * (2*mLi*kB*Tabs./...
    (TTilde_vec(length(TTilde_vec)-length(select_mu_TTilde_vec)+1:end)*hbar^2)).^(3/2);

close all

figure(1)
plot(z_TTilde_vec,n_vec_TTilde)

save('UnitaryFakeProfile_T0_1');