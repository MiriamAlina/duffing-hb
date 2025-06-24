%% Duffing Oscillator: Harmonic Balance Analysis
% Based on: Krack & Gross, 
%           Harmonic Balance for Nonlinear Vibration Problems, 
%           Chapter 5, Exercise 1
% Author:   Miriam Goldack
% Date:     2025-06-23

addpath(genpath('NLvib-NLvib-Basic'));


%% System Parameters (Section 5.2a)
mu    = 1;
kappa = 1;
zeta  = 0.05;
gamma = 0.1;
P     = 0.18;
H     = 7;     % Number of harmonics 
N     = 2^6;   % Number of time samples


%% Numerical Harmonic Balance Solution (Section 5.2b)
Om_s = 0.5;    % Start frequency
Om_e = 1.6;    % End frequency

% Initial guess (based on linearized system)
Q0 = (-Om_s^2 * mu + 1i * Om_s * zeta + kappa) \ P; 
x0 = [0; real(Q0); -imag(Q0); zeros(2*(H-1), 1)];

ds   = 0.01;                   % Step size
Sopt = struct('jac', 'none');  % Solver options

X = solve_and_continue(x0, @(X) HB_residual_Duffing(X, mu, zeta, kappa, ...
    gamma, P, H, N), Om_s, Om_e, ds, Sopt);
Om = X(end, :);


%% Amplitude of Fundamental Harmonic (Section 5.2c)
a = sqrt(X(2, :).^2 + X(3, :).^2);


%% Analytical Single-Term HB Approximation (Section 5.3)
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


%% Frequency Response Plot (Figure Left)
figure('Position', [10 10 1200 300]);
subplot(1,2,1); hold on;

plot(Om_ana(valid_ana, 1), a_ana(valid_ana), 'g-', ...
    'DisplayName', 'analytical HB, H=1');
plot(Om_ana(valid_ana, 2), a_ana(valid_ana), 'g-', ...
    'HandleVisibility', 'off');
plot(Om, a, 'k--', 'DisplayName', ['numerical HB, H=' num2str(H)]);

xlabel('excitation frequency $\Omega$', 'Interpreter', 'latex');
ylabel('response amplitude $a$', 'Interpreter', 'latex');
xlim([Om_s Om_e]);
legend('Location', 'northwest', 'Interpreter', 'latex');
title('Amplitude-frequency curve', 'Interpreter', 'latex')

% Inset Zoom
ax_inset = axes('Position', [0.38 0.45 0.07 0.3]); 
box on; hold on;
plot(Om_ana(valid_ana, :), a_ana(valid_ana), 'g-');
plot(Om, a, 'k--');
xlim([1.25 1.275]);
ylim([2.6 2.9]);


%% Time Domain Comparison at Resonance (Section 5.4)
[~, i_res] = max(a);
Xres       = X(:, i_res);
Om_res     = Om(i_res);

% Harmonic Coefficients
h     = 1:H;
Q_cos = Xres(2:2:end-1);    
Q_sin = Xres(3:2:end-1);   
Q_ce  = [flipud(Q_cos + 1i * Q_sin)/2; Xres(1); (Q_cos - 1i * Q_sin)/2];

% Response over one period
t = linspace(0, 2*pi / Om_res, N);  
E = exp(1i * 2*pi/N * (0:N-1)' * (-H:H));  
q = real(E * Q_ce);

% Derivatives
dt  = t(2) - t(1);
dq  = gradient(q, dt);
ddq = gradient(dq, dt);

% Analytical Time Response at Resonance from HB (H=1)
a_valid          = a_ana(valid_ana);
Om_valid         = Om_ana(valid_ana, :);
[~, idx_ana_res] = max(a_valid);
a_res_ana        = a_valid(idx_ana_res);
Om_res_ana       = Om_valid(idx_ana_res, 1);

a1 = 0;          % Cosine coefficient 
b1 = a_res_ana;  % Sine coefficient

q_ana   = a1 * cos(Om_res_ana * t) + b1 * sin(Om_res_ana * t);
dq_ana  = -a1 * Om_res_ana * sin(Om_res_ana * t) + ...
    b1 * Om_res_ana * cos(Om_res_ana * t);
ddq_ana = -a1 * Om_res_ana^2 * cos(Om_res_ana * t) - ...
    b1 * Om_res_ana^2 * sin(Om_res_ana * t);


%% Time Signals Plot (Figure Right)
subplot(1,2,2); hold on;
plot(t, q_ana,   'g-');
plot(t, q,       'k--');
plot(t, dq_ana,  'g-');
plot(t, dq,      'k--');
plot(t, ddq_ana, 'g-');
plot(t, ddq,     'k--');

xlabel('$$\Omega\, t$$', 'Interpreter', 'latex');
title(['Time history of $q, \dot{q}$ and $\ddot{q}$ at the point of ' ...
    'maximum amplitude'], 'Interpreter', 'latex')

% Annotations
annotation('textarrow', [0.70 0.675], [0.80 0.73], 'String', ...
    '$$q$$', 'Interpreter', 'latex');
annotation('textarrow', [0.62 0.600], [0.82 0.75], 'String', ...
    '$$\dot{q}$$', 'Interpreter', 'latex');
annotation('textarrow', [0.87 0.850], [0.87 0.80], 'String', ...
    '$$\ddot{q}$$', 'Interpreter', 'latex');

