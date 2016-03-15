[ KappaTilde_ideal, PTilde_ideal, TTilde_ideal, CV_NkB_ideal, beta_mu_vec_ideal ,Z_vec_ideal ] ...
    = IdealFermiEOS(10^7);

BetaMu_vec_ideal = beta_mu_vec_ideal(1) + ...
    cumtrapz((1./(KappaTilde_ideal.*(TTilde_ideal).^2)),TTilde_ideal);

beta_mu = zeros(1,length(TTilde_ideal));
beta_mu(1) = beta_mu_vec_ideal(1);
core = (1./(KappaTilde_ideal.*(TTilde_ideal).^2));
DT = diff(TTilde_ideal);
    for i=2:length(TTilde_ideal)-1
    beta_mu(length(TTilde_ideal)-(i-1)) = -  DT(i:end)*core(i+1:end)' + beta_mu(1);
    end

figure(2)
plot(TTilde_ideal,(beta_mu_vec_ideal-BetaMu_vec_ideal))
ylim([-0.5 0.5])

figure(3)
plot(TTilde_ideal,beta_mu_vec_ideal-beta_mu)
ylim([-0.01 0.01])

figure(4)
plot(TTilde_ideal,beta_mu)