%========================================================================
% DESCRIPTION: 
% The geometrically nonlinear, discretized Timoshenko beam under harmonic
% forcing at its tip is considered. The model is derived from the 
% publically available SSMTool (https://github.com/jain-shobhit/SSMTool),
% and consist of 4 elements with a total of 21 (nonlinear) degrees of 
% freedom. Material properties and forcing pattern is set according to the
% numerical example in Ponsioen et al., JSV, 2020 [1].
% 
% NOTE: Unit system may be inconsistent, length-force-weight: [mm-mN-kg]
% is used in the SSMTool. In this system elastic modulus 'E' should be: 
% E=90GPa==90.000 [N/(mm^2)]==90e6 [mN/(mm^2)]; in FEM_Timoshenko.m 
% E=90.000 [mN/(mm^2)] is used. In this code, we stick to the numerical 
% values from the paper.
%
% In [1], NLvib was used and found to struggle: It took 3 days to obtain 
% 87,340 solution points which span only a small portion of the frequency 
% response curve shown in Figure 6.
% 
% In this example, we show that the complete curve can be computed with ca.
% 100 solution points for an appropriate choice of parameters, which should
% run through in less than 2 minutes on a standard computer. We approached
% the authors of [1] to obtain the NLvib settings they used to reproduce the
% poor performance. Apparently, the main problem was that they did not use
% appropriate values for the DSCALE parameter and did not make use of the
% dynamicDSCALE option. The sensitivity to those parameters is well-explained
% in the manual and the book. Details are discussed in [2].
%
% We acknowledge that it may require experience to set parameters of 
% numerical tools properly. We are always happy to try to help if you 
% encounter any problems, which you cannot solve with the manual or the book.
% 
% The HB-approximation is verified by numerical time integration using
% ode15s.
% 
% [1] https://doi.org/10.1016/j.jsv.2020.115640
% [2] https://doi.org/10.1016/j.jsv.2024.118808
%========================================================================
% This file is part of NLvib.
% 
% If you use NLvib, please refer to the book:
%   M. Krack, J. Gross: Harmonic Balance for Nonlinear Vibration
%   Problems. Springer, 2019. https://doi.org/10.1007/978-3-030-14023-6.
% 
% COPYRIGHT AND LICENSING: 
% NLvib Version 1.3 Copyright (C) 2020  Malte Krack  
%										(malte.krack@ila.uni-stuttgart.de) 
%                     					Johann Gross 
%										(johann.gross@ila.uni-stuttgart.de)
%                     					University of Stuttgart
% This program comes with ABSOLUTELY NO WARRANTY. 
% NLvib is free software, you can redistribute and/or modify it under the
% GNU General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% For details on license and warranty, see http://www.gnu.org/licenses
% or gpl-3.0.txt.
%========================================================================
clearvars;
close all;
clc;

addpath('../../SRC');
addpath('../../SRC/MechanicalSystems');
%% Build model from [1]
% Geometrically nonlinear Timoshenko beam
loadModel = true; % set to false if you vary parameters
if ~loadModel
    addpath('EXT/SSMTool');
    nElements = 4;
    isViscoelastic = 0;
    [M,D,K,~,outdof,p,E] = ...
        build_model(nElements,isViscoelastic);
    save('Timoshenko.mat','M','D','K','outdof','p','E');
else
    load('Timoshenko.mat','M','D','K','outdof','p','E');
end

% Set forcing level according to [1] to 2.4kN
Fex1  = zeros(size(M,1),1); Fex1(end-1) = 2.4e6;

% Set response locator vector
Tresp = zeros(size(M,1),1)'; Tresp(outdof) = 1;

% Define oscillator as system with polynomial stiffness nonlinearities
timoshenkoBeam = ...
    System_with_PolynomialStiffnessNonlinearity(M,D,K,p,E,Fex1);

% Number of degrees of freedom
n = timoshenkoBeam.n;

% Frequency of first linear mode
[phi1,om12] = eigs(K,M,1,'sm'); % = 3.05 (plausible, cf. Fig. 6)
om1 = sqrt(om12);

%% Compute frequency response using harmonic balance

% Analysis parameters
analysis = 'FRF';
H = 10;      % harmonic order
N = 4*H+1;  % number of samples sufficient for max. cubic-degree polynomials, cf. Appendix A in https://doi.org/10.1016/j.ymssp.2019.106503
Om_s = 2.5; % start frequency
Om_e = 3.9; % end frequency

% Initial guess (solution of underlying linear system)
Q1 = (-Om_s^2*M + 1i*Om_s*D + K)\timoshenkoBeam.Fex1;
y0 = zeros((2*H+1)*length(Q1),1);
y0(length(Q1)+(1:2*length(Q1))) = [real(Q1);-imag(Q1)];
qscl = max(abs(y0));

% Solve and continue w.r.t. Om
ds = .05;
Sopt = struct('Dscale',[qscl*ones(size(y0));Om_s],'dynamicDscale',1);
tmp = tic;
X_HB = solve_and_continue(y0,...
    @(X) HB_residual(X,timoshenkoBeam,H,N,analysis),...
    Om_s,Om_e,ds,Sopt);
tFRC = toc(tmp);

% Interpret solver output
Om_HB = X_HB(end,:);
% convert Fourier coefficients to complex exponential notation
I0          = 1:n; ID = n+(1:H*n);
IC          = n+repmat(1:n,1,H)+n*kron(0:2:2*(H-1),ones(1,n)); IS = IC+n;
Q           = zeros(n*(H+1),size(X_HB,2));
Q(I0,:)     = X_HB(I0,:);
Q(ID,:)     = X_HB(IC,:)-1i*X_HB(IS,:);
% Redefine amplitude as maximum amplitude at response DOF
tau         = linspace(0,2*pi,2^7)';
H_iDFT      = exp(1i*tau*(0:H));
Qresp       = kron(eye(H+1),Tresp)*Q;
qFRC_max    = zeros(size(Qresp,2),1);
for nn=1:length(qFRC_max)
    qFRC_max(nn) = max((real(H_iDFT*Qresp(:,nn))));
end

figure(1);
plot(Om_HB,qFRC_max,'b-','LineWidth',2,'DisplayName',...
    sprintf('H=%d, computed in %.2f min',H,round(tFRC/60,2)));
hold on;
set(gca,'xlim',[2.5 3.9],'ylim',[0 .8]);
xlabel('\Omega');
ylabel('max(q_{resp})');
legend()

%% Compute frequency response using time step integration

% Exponents and coefficients of the geom. nonlinear terms
nlElts  = timoshenkoBeam.nonlinear_elements{1};
pp      = nlElts.exponents;
Et      = transpose(nlElts.coefficients);
nz      = size(Et,2);
% define function handle for nonlinear terms evaluation
fnl     = @(q) (Et*reshape(prod(kron(q',ones(size(pp,1),1)).^pp,2),nz,1))';

% Define right hand side of ODE
ODE = @(t,y,Om) [y(n+1:2*n);M\( Fex1*cos(Om*t) - K*y(1:n) - D*y(n+1:2*n) ...
    -transpose(fnl(y(1:n))) )];

% Exclude overhanging branch (presuming Om_s<Om_e here)
valid = find(diff(Om_HB)>0);
% Select a subset of solution points for validation (here every 4th pt.)
valid = valid(1:4:end);
Om_TI       = Om_HB(valid);
Q_TI        = Q(:,valid);
% Indicate selected points
plot(Om_HB(valid),qFRC_max(valid),'bo','MarkerSize',6, ...
    'DisplayName','selected')

% Loop over selected points
qTI_max     = zeros(length(Om_TI),1);
for iom =1:length(Om_TI)
    fprintf('Time integration for point %d/%d is running.\n',...
        iom,length(Om_TI))

    % Select initial conditions based on HB prediction
    omTmp  = Om_TI(iom);        
    q0     =  real(kron(exp(1i*tau(1)*(0:H)),eye(n))*Q_TI(:,iom));
    dot_q0 =  real(kron(exp(1i*tau(1)*(0:H)).*(1i*omTmp*(0:H)),eye(n))*...
        Q_TI(:,iom));
    
    % Simulate over 5 excitation periods (5 may be a bit short in general;
    % here it seems OK since damping is relatively high and the HB
    % approximation is pretty accurate)
    nPer   = 5;
    t_sim  = nPer*2*pi/omTmp;
    sol    = ode15s(@(t,y) ODE(t,y,omTmp),[0 t_sim],[q0;dot_q0]);

    % Evaluate response level in last period, use 50 samples per period for
    % post-processing
    nSamp = 50;
    n_steps = nSamp*round(t_sim*(omTmp/2/pi));
    t = linspace(0,t_sim,n_steps);
    y = deval(sol,t);        
    qResp = y(outdof,:);       
    qTI_max(iom) = max(qResp(end-nSamp:end));

%     % visual check whether steady state is reached
%     figure(2);
%     plot(t,qResp,'r'); hold on;
%     plot(t(end-nSamp:end),qResp(end-nSamp:end),'b','LineWidth',2);hold off
%     pause(.1)

end

figure(1); hold on
plot(Om_HB(valid),qTI_max,'rx','MarkerSize',6, ...
    'DisplayName','time integration')
legend('Location','best')
disp('finished')
