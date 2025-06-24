%% Duffing Oscillator: Harmonic Balance Analysis
% Based on: Krack & Gross, 
%           Harmonic Balance for Nonlinear Vibration Problems, 
%           Chapter 1, Introduction, Fig. 1.2
% Author:   Miriam Goldack
% Date:     2025-06-23

figure('Position', [10 10 600 400]);hold on;
Om_s = 0.2;
Om_e = 2.5;
Sopt = struct('jac', 'none');

colorList = lines(4); 
P_values = [0.02, 0.05, 0.1, 0.2];

for i = 1:length(P_values)
    P_level = P_values(i);
    system_settings = struct('mu', 1, ...
        'kappa', 1, ...
        'zeta', 0.05, ...
        'gamma', 1, ...
        'P', P_level, ...
        'H', 7, ...
        'N', 2^6, ...
        'Om_s', 0.2, ...
        'Om_e', 2.5);
    analysis_result = HB_analysis(system_settings, Sopt);
    Om_ana_valid = analysis_result.Om_ana_valid;
    a_ana_valid = analysis_result.a_ana_valid;
    Om = analysis_result.Om;
    a = analysis_result.a;
    
    color = colorList(i, :);
    plot(Om_ana_valid(:,1), a_ana_valid, '-', 'Color', color, ...
        'DisplayName', ['H=1, P=' num2str(P_level)]);
    plot(Om_ana_valid(:,2), a_ana_valid, '-', 'Color', color, ...
        'HandleVisibility', 'off');
    plot(Om, a, '.', 'Color', color, ...
        'DisplayName', ['H=' num2str(H) ', P=' num2str(P_level)]);    
end

xlabel('excitation frequency $\Omega$', 'Interpreter', 'latex');
ylabel('response amplitude $a$', 'Interpreter', 'latex');
xlim([Om_s Om_e]);
legend('Location', 'northwest', 'Interpreter', 'latex');
title(['Amplitude-frequency curves for different levels of ' ...
    'excitation $P$'], 'Interpreter', 'latex')