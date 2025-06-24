%% Duffing Oscillator: Harmonic Balance Analysis
% Based on: Krack & Gross, 
%           Harmonic Balance for Nonlinear Vibration Problems, 
%           Chapter 5, Exercise 1
% Author:   Miriam Goldack
% Date:     2025-06-23


%% Harmonic Balance Residual for the Duffing Oscillator (Section 5.1)
function R = HB_residual_Duffing(X, mu, zeta, kappa, gamma,P, H, N)
    % INPUTS:
    %     X      - State vector [q_mean; cos_harmonics; sin_harmonics; Omega]
    %     mu     - Mass of the oscillator
    %     zeta   - Damping coefficient
    %     kappa  - Linear stiffness
    %     gamma  - Cubic (Duffing) stiffness coefficient
    %     P      - Amplitude of external cosine forcing
    %     H      - Number of harmonics used in the approximation
    %     N      - Number of time samples for Fourier transforms
    %
    % OUTPUT:
    %     R      - Residual vector in sine/cosine formulation (real-valued)
    
    % Conversion of sine-cosine to complex-exponential representation
    Q_cos = X(2:2:end-1);  % Cosine terms of harmonics
    Q_sin = X(3:2:end-1);  % Sine terms of harmonics
    Q_ce = [flipud(Q_cos + 1i * Q_sin)/2; ...  % Negative harmonics
            X(1); ...                          % Mean value (DC)
            (Q_cos - 1i * Q_sin)/2];           % Positive harmonics
    
    % Excitation frequency
    Om = X(end);  
    
    % Forcing in complex exponential form (cosine forcing)
    Fex_ce = [zeros(H-1, 1); P/2; 0; P/2; zeros(H-1, 1)];
    
    % Inverse DFT matrix (time domain reconstruction)
    E_NH = exp(1i * 2*pi/N * (0:N-1)' * (-H:H));
    
    % Time-domain response q(t)
    q = real(E_NH * Q_ce);
    
    % Nonlinear force in time domain (Duffing nonlinearity)
    fnl = gamma * q.^3;
    
    % Fourier coefficients of nonlinear force (DFT)
    Fnl_ce = E_NH' / N * fnl;
    
    % Dynamic equilibrium in frequency domain
    omega_h = (-H:H)' * Om;  % Harmonic frequencies
    R_ce = (-mu * omega_h.^2 + 1i * omega_h * zeta + kappa) .* Q_ce + ...
        Fnl_ce - Fex_ce;
    
    % Conversion of complex-exponential to sine-cosine representation
    R = [real(R_ce(H+1)); ...      % Mean value
         real(R_ce(H+2:end)); ...  % Cosine coefficients
        -imag(R_ce(H+2:end))];     % Sine coefficients
end