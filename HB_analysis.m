%% Duffing Oscillator: Harmonic Balance Analysis
% Based on: Krack & Gross, 
%           Harmonic Balance for Nonlinear Vibration Problems,
%           Chapter 5, Exercise 1
% Author:   Miriam Goldack
% Date:     2025-06-23


function analysis_result = HB_analysis(system_settings, Sopt)
    mu = system_settings.mu;
    kappa = system_settings.kappa;
    zeta = system_settings.zeta;
    gamma = system_settings.gamma;
    P = system_settings.P;
    H = system_settings.H;
    N = system_settings.N;
    Om_s = system_settings.Om_s;
    Om_e = system_settings.Om_e;

    % Initial guess (based on linearized system)
    Q0 = (-Om_s^2 * mu + 1i * Om_s * zeta + kappa) \ P; 
    x0 = [0; real(Q0); -imag(Q0); zeros(2*(H-1), 1)];
    
    ds   = 0.01;                    
    
    X = solve_and_continue(x0, @(X) HB_residual_Duffing(X, mu, zeta, kappa, ...
        gamma, P, H, N), Om_s, Om_e, ds, Sopt);
    Om = X(end, :);
    
    a = sqrt(X(2, :).^2 + X(3, :).^2);
    
    a_min  = 0.1;
    a_max  = 3;
    a_ana  = linspace(a_min, a_max, 50);
    Om_ana = zeros(length(a_ana), 2);
    
    for i = 1:length(a_ana)
        a_sqr      = a_ana(i)^2;
        inner_sqrt = sqrt(P^2 / a_sqr + zeta^4 / 4 - zeta^2 - ...
            (3 * zeta^2 * gamma * a_sqr) / 4);
        base       = 1 - zeta^2 / 2 + (3 * gamma * a_sqr) / 4;
    
        Om_ana(i,1) = sqrt(base + inner_sqrt);
        Om_ana(i,2) = sqrt(base - inner_sqrt);
    end
    
    valid_ana = imag(Om_ana(:, 1)) == 0 & imag(Om_ana(:, 2)) == 0;
    Om_ana_valid = Om_ana(valid_ana, :);
    a_ana_valid = a_ana(valid_ana);

    analysis_result = struct('Om_ana_valid', Om_ana_valid, ...
        'a_ana_valid', a_ana_valid, ...
        'Om', Om, ...
        'a', a);
end