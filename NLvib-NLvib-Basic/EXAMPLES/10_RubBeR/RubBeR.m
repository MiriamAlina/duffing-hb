%========================================================================
% DESCRIPTION: 
% Investigation of the dynamics of a cantilevered Euler-Bernoulli beam
% with Jenkins element, similar (!) to the one experimentally studied in 
% Scheel et al., JSV, 2020. doi: 10.1016/j.jsv.2020.115580.
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
%% Define system

% Properties of the beam
len = .710;             % length
height = .05;           % height in the bending direction
thickness = .06;        % thickness in the third dimension
E = 210e9*0.6725;       % Young's modulus
rho = 7850;             % density
BCs = 'clamped-free';   % constraints

% Set up one-dimensional finite element model of an Euler-Bernoulli beam:
% Discretization of RubBeR free length into 14 elements of 50 mm each.
n_nodes = 15;           % number of equidistant nodes along length
beam = FE_EulerBernoulliBeam(len,height,thickness,E,rho,...
    BCs,n_nodes);

% Specify elastic dry friction (Jenkins) element
kt = 1e8;               % tangential stiffness
muN = .5*53*4;          % limit friction force
inode = 6;              % ~2nd pocket
dir = 'trans';          % transverse direction
add_nonlinear_attachment(beam,inode,dir,'elasticDryFriction',...
    'stiffness',kt,'friction_limit_force',muN,'ishysteretic',1);

% Vectors recovering deflection at tip and nonlinearity location
T_nl = beam.nonlinear_elements{1}.force_direction';
%% Modal analysis the linearized system

% Modes for free sliding contact
[PHI_free,OM2] = eig(beam.K,beam.M);
om_free = sqrt(diag(OM2));
% Sorting
[om_free,ind] = sort(om_free); PHI_free = PHI_free(:,ind);

% Modes for fixed contact
inl = find(T_nl); w = zeros(length(beam.K),1); w(inl) = 1;
Kt = kt*(w*w');
[PHI_stick,OM2] = eig(beam.K+Kt,beam.M);
om_stick = sqrt(diag(OM2));
% Sorting
[om_stick,ind] = sort(om_stick); PHI_stick = PHI_stick(:,ind);

%% Nonlinear modal analysis using harmonic balance
analysis = 'NMA';

% Analysis parameters
imod = 1;           % mode to be analyzed
H = 7;              % harmonic order
N = 2^8;           % number of time samples per period
log10a_s = -5;      % start vibration level (log10 of modal mass)
log10a_e = -2.5;    % end vibration level (log10 of modal mass)
inorm = beam.n-1;   % coordinate for phase normalization

% Initial guess vector x0 = [Psi;om;del], where del is the modal
% damping ratio, estimate from underlying linear system (sticking)
model = beam; n = size(PHI_stick,1); 
om = om_stick(imod); phi = PHI_stick(:,imod);
Psi = zeros((2*H+1)*n,1);
Psi(n+(1:n)) = phi;
x0 = [Psi;om;0];
psiscl = max(abs((Psi)));

% Solve and continue w.r.t. Om
ds = .02;
Sopt = struct('flag',1,'stepadapt',1,'dynamicDscale',1,...
    'Dscale',[1e-0*psiscl*ones(size(Psi));...
    (om_stick(imod)+om_free(imod))/2;1e-1;1e0]);
fscl = mean(abs(model.K*phi));
X_NM = solve_and_continue(x0,...
    @(X) HB_residual(X,model,H,N,analysis,inorm,fscl),...
    log10a_s,log10a_e,ds,Sopt);

% Interpret solver output
Psi_NM = X_NM(1:end-3,:);
om_NM = X_NM(end-2,:);
del_NM = X_NM(end-1,:);
log10a_NM = X_NM(end,:);
a_NM = 10.^log10a_NM;
Q_NM = Psi_NM.*repmat(a_NM,size(Psi_NM,1),1);
Psi1_NM = (Psi_NM(n+(1:n),:)-1i*Psi_NM(2*n+(1:n),:));
%% Amplitude-depndent modal frequency and damping ratio

% Modal frequency vs. amplitude
figure(1); hold on;
plot(a_NM,om_NM/om_stick(imod),'k-');
set(gca,'xscale','log','xlim',[min(a_NM) max(a_NM)]);
xlabel('modal amplitude'); ylabel('\omega/\omega_{lin}');

% Modal damping ratio vs. amplitude
figure(2); hold on;
plot(a_NM,del_NM*1e2,'k-');
set(gca,'xscale','log','xlim',[min(a_NM) max(a_NM)]);
xlabel('modal amplitude'); ylabel('modal damping ratio in %');
%% Mode complexity at maximum damping point
figure(3);
[~,imaxD] = max(del_NM);
tmp = Psi1_NM(1:2:end,imaxD);       % select transferse displacements
tmp = tmp*abs(tmp(end))/tmp(end);   % phase normalize to tip
tmp = tmp./abs(tmp);                % length normalize to 1
compass(tmp);
xlabel('real part of norm. modal deflection shape');
ylabel('imaginary part of norm. modal deflection shape');
%% Hrmonic and modal contributions
E = zeros(n,H,size(om_NM,2));
for h=1:H
    Qh = Q_NM(n*(2*h-1)+(1:n),:) - ...
        1i*Q_NM(n*(2*h-1)+n+(1:n),:);
    Qhmod = PHI_stick'*beam.M*Qh;
    E(:,h,:) = 1/4*( (h*repmat(om_NM,n,1)).^2 + ...
        repmat(om_stick(:),1,size(om_NM,2)).^2 ) .* abs(Qhmod).^2;
end
% Identify most important contributions
[row,col] = find(sum(E./repmat(E(1,1,:),n,H,1)>1e-2,3)>=1);
ref = E(1,1,:);

figure; hold on;
legstr = cell(length(row),1);
for i=2:length(row)
    plot(a_NM,squeeze(E(row(i),col(i),:)./ref));
    legstr{i} = ['m=' num2str(row(i)) ',h=' num2str(col(i))];
end
legend(legstr{:});
set(gca,'xscale','log','yscale','log','ylim',[1*10^-2 10],...
    'xlim',[min(a_NM) max(a_NM)]);
xlabel('modal amplitude'); ylabel('E(m,h)/E(1,1)');