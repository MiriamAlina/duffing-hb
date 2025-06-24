%% Duffing Oscillator Revisited: Harmonic Balance Analysis
% Based on: Krack & Gross, 
%           Harmonic Balance for Nonlinear Vibration Problems, 
%           Chapter 5, Homework B, 5.17 Going Around Turning Points
% Author:   Miriam Goldack
% Date:     2025-06-23


%% (a), (b), (c) 
system_settings = struct('mu', 1, ...
    'kappa', 1, ...
    'zeta', 0.01, ...
    'gamma', 0.1, ...
    'P', 1, ...
    'H', 7, ...
    'N', 2^6, ...
    'Om_s', 0.8, ...
    'Om_e', 1.2);
Sopt = struct('jac', 'none');
analysis_result = HB_analysis(system_settings, Sopt);
Om_ana_valid = analysis_result.Om_ana_valid;
a_ana_valid = analysis_result.a_ana_valid;
Om = analysis_result.Om;
a = analysis_result.a;

figure('Position', [10 10 1600 300]);
subplot(1,5,1); hold on;

plot(Om_ana_valid(:,1), a_ana_valid, 'g-', ...
    'DisplayName', 'analytical HB, H=1');
plot(Om_ana_valid(:,2), a_ana_valid, 'g-', ...
    'HandleVisibility', 'off');
plot(Om, a, 'k--', 'DisplayName', ['numerical HB, H=' num2str(H)]);

xlabel('excitation frequency $\Omega$', 'Interpreter', 'latex');
ylabel('response amplitude $a$', 'Interpreter', 'latex');
legend('Location', 'northwest', 'Interpreter', 'latex');
title('5.17 (c) Little damping, strong excitation', ['continuation ' ...
    'only for stable small frequencies'], 'Interpreter', 'latex')


%% (d)
system_settings.Om_s = 1.2;    % Start frequency
system_settings.Om_e = 0.8;    % End frequency
disp(system_settings)

% analysis_result = HB_analysis(system_settings, Sopt);

disp('--------------------')
disp('5.17 (d)')
disp('After swapping the frequency limits, the output of NVlib is this:')
disp('"Provided initial guess is not in the basin of attraction."')


%% (e)
Sopt.reversaltolerance = inf; 

% analysis_result = HB_analysis(system_settings, Sopt);

disp('--------------------')
disp('5.17 (e)')
disp('Setting reversaltolerance to inf still gives the same error:')
disp('"Provided initial guess is not in the basin of attraction."')


%% (f)
system_settings.Om_s = 0.8;    
system_settings.Om_e = 1.8;  

analysis_result = HB_analysis(system_settings, Sopt);
Om_ana_valid = analysis_result.Om_ana_valid;
a_ana_valid = analysis_result.a_ana_valid;
Om = analysis_result.Om;
a = analysis_result.a;

subplot(1,5,2); hold on;

plot(Om_ana_valid(:, 1), a_ana_valid, 'g-', ...
   'DisplayName', 'analytical HB, H=1');
plot(Om_ana_valid(:, 2), a_ana_valid, 'g-', ...
   'HandleVisibility', 'off');
plot(Om, a, 'k--', 'DisplayName', ['numerical HB, H=' num2str(H)]);

xlabel('excitation frequency $\Omega$', 'Interpreter', 'latex');
ylabel('response amplitude $a$', 'Interpreter', 'latex');
legend('Location', 'northwest', 'Interpreter', 'latex');
title('5.17 (f) Turning point bypass', ['only following the left ' ...
    'branch'], 'Interpreter', 'latex')


%% (g)
Sopt.flag = 0;

analysis_result = HB_analysis(system_settings, Sopt);
Om_ana_valid = analysis_result.Om_ana_valid;
a_ana_valid = analysis_result.a_ana_valid;
Om = analysis_result.Om;
a = analysis_result.a;

subplot(1,5,3); hold on;

plot(Om_ana_valid(:, 1), a_ana_valid, 'g-', ...
    'DisplayName', 'analytical HB, H=1');
plot(Om_ana_valid(:, 2), a_ana_valid, 'g-', ...
    'HandleVisibility', 'off');
plot(Om, a, 'k--', 'DisplayName', ['numerical HB, H=' num2str(H)]);

xlabel('excitation frequency $\Omega$', 'Interpreter', 'latex');
ylabel('response amplitude $a$', 'Interpreter', 'latex');
legend('Location', 'northwest', 'Interpreter', 'latex');
title('5.17 (g) Continuation off', 'remaining on the same branch', ...
    'Interpreter', 'latex')

 
%% (h) 
system_settings.Om_s = 1.8;  
system_settings.Om_e = 0.8;    

analysis_result = HB_analysis(system_settings, Sopt);
Om_ana_valid = analysis_result.Om_ana_valid;
a_ana_valid = analysis_result.a_ana_valid;
Om = analysis_result.Om;
a = analysis_result.a;

subplot(1,5,4); hold on;

plot(Om_ana_valid(:,1), a_ana_valid, 'g-', ...
    'DisplayName', 'analytical HB, H=1');
plot(Om_ana_valid(:,2), a_ana_valid, 'g-', ...
    'HandleVisibility', 'off');
plot(Om, a, 'k--', 'DisplayName', ['numerical HB, H=' num2str(H)]);

xlabel('excitation frequency $\Omega$', 'Interpreter', 'latex');
ylabel('response amplitude $a$', 'Interpreter', 'latex');
legend('Location', 'northwest', 'Interpreter', 'latex');
title('5.17 (h) Starting point at higher frequency', ['convergence to ' ...
    'right branch'], 'Interpreter', 'latex')
 
 
%% (i)
system_settings.P = 0.05;

analysis_result = HB_analysis(system_settings, Sopt);
Om_ana_valid = analysis_result.Om_ana_valid;
a_ana_valid = analysis_result.a_ana_valid;
Om = analysis_result.Om;
a = analysis_result.a;

subplot(1,5,5); hold on;

plot(Om_ana_valid(:,1), a_ana_valid, 'g-', ...
    'DisplayName', 'analytical HB, H=1');
plot(Om_ana_valid(:,2), a_ana_valid, 'g-', ...
    'HandleVisibility', 'off');
plot(Om, a, 'k--', 'DisplayName', ['numerical HB, H=' num2str(H)]);

xlabel('excitation frequency $\Omega$', 'Interpreter', 'latex');
ylabel('response amplitude $a$', 'Interpreter', 'latex');
legend('Location', 'northwest', 'Interpreter', 'latex');
title('5.17 (i) Small excitation', 'insufficient numerical resolution', ...
    'Interpreter', 'latex')