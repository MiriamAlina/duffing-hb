%========================================================================
% DESCRIPTION: 
% A linear clamped-free Euler-Bernoulli beam is considered with a cubic 
% spring attached to its free end. The beam is discretized by standard 
% beam finite elements. The number of degrees of freedom can be varied by
% adjusting the number of finite elements. The problem is considered in
% [1], where it was concluded that Harmonic Balance leads to intractable
% comptuation effort.
% 
% An important aspect of the problem is that the nonlinear element, i.e., 
% the cubic spring, is associated with a single degree of freedom. In other 
% words, the nonlinear force depends on and acts on only a single
% coordinate. For brevity, we refer to this coordinate as nonlinear 
% coordinate, and the remaining ones as linear coordinates. If the number 
% of nonlinear coordinates is much smaller than that of the linear 
% coordinates, we call the nonlinear terms sparse. The sparsity is 
% inherited by the Harmonic Balance equations. One can exploit this 
% sparsity by an exact condensation. Essentially, one eliminates the 
% Fourier coefficients of the linear coordinates within the algebraic 
% equation system.
% 
% We have implemented the particular variant of the exact condensation 
% described in the HB book (Section 4.3), and is part of the Homework 
% Problem F (Section 5). With this, the computation effort does not seem
% intractable. Details are discussed in [2].
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
% NLvib Version 1.3 Copyright (C) 2024  Malte Krack  
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
%% Set number of degrees of freedom 'n' and resulting number of nodes
n = 30; % [10||30|50|70|100|200|500|1000|10000]
n_nodes = n/2+1;
%% Define finite element model

% Properties of the beam
len = 2.7;              % length
height = 10e-3;         % height in the bending direction
thickness = 10e-3;      % thickness in the third dimension
E = 45e9;               % Young's modulus
rho = 1780;             % density
BCs = 'clamped-free';   % constraints

% Setup one-dimensional finite element model of an Euler-Bernoulli beam
beamFE = FE_EulerBernoulliBeam(len,height,thickness,E,rho,...
    BCs,n_nodes);

% Use sparse algebra
beamFE.M = sparse(beamFE.M);
beamFE.K = sparse(beamFE.K);
beamFE.Fex1 = sparse(beamFE.Fex1);

% Apply cubic element at free end in translational direction
inode = n_nodes;
dir = 'trans';
kappa = 4e6*1e-6;
add_nonlinear_attachment(beamFE,inode,dir,'cubicSpring',...
    'stiffness',kappa);

% Define damping
alpha = 1.25e-4;
beta   = 2.5e-4;
beamFE.D = alpha*beamFE.M+beta*beamFE.K;

% Apply forcing to free end of beam in translational direction, with
% magnitude 'fex'
inode = n_nodes;
dir = 'trans';
fex = 2e-3;
add_forcing(beamFE,inode,dir,fex);

%% Prepare dynamic condensation to the nonlinear part
disp('START preparing dynamic condensation');
tic;

% This mechanical system class is specific to the condensation
beam = MechanicalSystem_nlred(beamFE.M,beamFE.D,beamFE.K,...
    beamFE.nonlinear_elements,beamFE.Fex1);
n = length(beam.M);

% Compute linear modes (required for spectral decomposition of dynamic
% compliance matrix, which is used for condensation)
% FOR NUMERICAL STABILITY AT HIGH #DOFs, WE COMPUTE LOWEST-FREQUENCY MODES
% SEPARATELY
tic
if n>200
    Nu = 100;
else
    Nu = 1;
end
[PHI_,LAM_] = eigs(beam.K,beam.M,Nu,'smallestabs');
[om_,ind] = sort(sqrt(diag(LAM_)));
PHI_ = PHI_(:,ind);
PHI_ = PHI_./repmat(sqrt(diag(PHI_'*beam.M*PHI_))',size(PHI_,1),1);
[PHI__,LAM__] = eigs(beam.K,beam.M,n-Nu,'largestabs');
[om__,ind] = sort(sqrt(diag(LAM__)));
PHI__ = PHI__(:,ind);
PHI__ = PHI__./repmat(sqrt(diag(PHI__'*beam.M*PHI__))',size(PHI__,1),1);
PHI = [PHI_ PHI__];
om = [om_;om__];
disp(['Linear modes computed in ' num2str(toc) ' s.']);
disp(['Lowest natural freq. is ' num2str(om(1)) ' rad/s.']);

% Determine indices associated with nonlinearities
W = zeros(n,length(beam.nonlinear_elements));
for i=1:length(beam.nonlinear_elements)
    if ~beam.nonlinear_elements{i}.islocal
        error(['Reduction to nonlinear coordinates only reasonable ' ...
            'if global nonlinearities are absent.']);
    end
    W(:,i) = beam.nonlinear_elements{i}.force_direction;
end
iN = find(sum(W,2));

% Store this information
beam.iN = iN;
beam.PHI = PHI;
beam.mmod = ones(1,n);
beam.kmod = (om(:).').^2;
beam.dmod = alpha + beta*beam.kmod;
disp('END preparing dynamic condensation');
%% Harmonic Balance computation WITH CONDENSATION

% Analysis parameters
analysis = 'FRF';
H = 10;         % harmonic order
N = 4*H+1;      % number of time samples per period cf. Appendix A in https://doi.org/10.1016/j.ymssp.2019.106503
Om_s = 6.88;    % start frequency
Om_e = 7.12;    % end frequency

% Initial guess (from underlying linear system)
Q1red = beam.PHI(iN,:)*diag(1./(-Om_s^2*beam.mmod + ...
    1i*Om_s*beam.dmod + beam.kmod))*(beam.PHI'*beam.Fex1);
qscl = mean(abs(Q1red));
x0 = zeros((2*H+1)*size(Q1red,1),1);
x0(size(Q1red,1)+(1:2*size(Q1red,1))) = [real(Q1red);-imag(Q1red)];
fscl = 1e-3;

% Solve and continue w.r.t. Om
ds = .5;
Sopt = struct('Dscale',[qscl*ones(size(x0));Om_s],'dynamicDscale',1);
X = solve_and_continue(x0,...
    @(X) HB_residual_nlred(X,beam,H,N,analysis,[],fscl),...
    Om_s,Om_e,ds,Sopt);

% Interpret solver output
Om = X(end,:);
Qtip = [X(1,:);X(2:2:end-1,:)-1i*X(3:2:end-1,:)];

% Determine absolute maximum value reached during one period
tau = linspace(0,2*pi,2^7)';
Qtip_max = max((real(exp(1i*tau*(0:H))*Qtip)));

% Illustrate frequency response
figure; hold on;
plot(Om,Qtip_max,'g-');
set(gca,'xlim',sort([Om_s Om_e]),'ylim',[0 .2],...
    'xtick',6.9:.05:7.1,'ytick',0:.02:.2);
grid on;
xlabel('\Omega'); 
ylabel('x_{n-1}');

%% Harmonic Balance computation WITHOUT CONDENSATION

% Run only for reasonably low 'n'
if n<100
    % Initial guess (from underlying linear system)
    Q1 = beam.PHI*diag(1./(-Om_s^2*beam.mmod + ...
        1i*Om_s*beam.dmod + beam.kmod))*(beam.PHI'*beam.Fex1);
    qscl = mean(abs(Q1));
    x0 = zeros((2*H+1)*size(Q1,1),1);
    x0(size(Q1,1)+(1:2*size(Q1,1))) = [real(Q1);-imag(Q1)];

    % Solve and continue w.r.t. Om
    ds = .5;
    Sopt = struct('Dscale',[4*sqrt(n)*qscl*ones(size(x0));Om_s],...
        'dynamicDscale',1,'eps',1e-3);
    X_ref = solve_and_continue(x0,...
        @(X) HB_residual(X,beamFE,H,N,analysis,[],fscl),...
        Om_s,Om_e,ds,Sopt);

    % Interpret solver output
    Om_ref = X_ref(end,:);
    Qtip_ref = [X_ref(iN,:);...
        X_ref(iN+n*(1:2:2*H),:)-1i*X_ref(iN+n*(2:2:2*H),:)];

    % Determine absolute maximum value reached during one period
    Qtip_max_ref = max((real(exp(1i*tau*(0:H))*Qtip_ref)));

    plot(Om_ref,Qtip_max_ref,'k--x');
    legend('HB w/ condensation','HB w/o condensation','LOCATION','BEST');
else
    legend('HB w/ condensation','LOCATION','BEST');
end
