%========================================================================
% DESCRIPTION: 
% Variant of HB_residual which implements an exact condensation to the
% nonlinear part. Global nonlinearities are presumed to be absent. The
% condensation relies on a spectral decomposition of the dynamic comliance
% matrix, see HB book, Section 4.3.
%========================================================================
function [R,dR,Q] = ...
    HB_residual_nlred(X,system,H,N,analysis_type,varargin)
%% Handle input variables depending on the modus

% 1 x n row vectors containing modal mass, damping and stiffness
mmod = system.mmod;
dmod = system.dmod;
kmod = system.kmod;
% PHI has the (mass-normalized) modal deflection shapes of the underlying 
% linear conservative system as columns (n x n matrix)
PHI = system.PHI;
% Indices of nonlinear coordinate(s)
iN = system.iN;

% number of degrees of freedom
n = length(mmod);
nN = length(iN);

% Conversion of real-valued to complex-valued harmonics of generalized
% coordinates q
I0 = 1:nN; ID = nN+(1:H*nN);
IC = nN+repmat(1:nN,1,H)+nN*kron(0:2:2*(H-1),ones(1,nN)); IS = IC+nN;
dX = eye(length(X));
QN = zeros(nN*(H+1),1);  dQN = zeros(size(QN,1),size(dX,2));
QN(I0) = X(I0);          dQN(I0,:) = dX(I0,:);
QN(ID) = X(IC)-1i*X(IS); dQN(ID,:) = dX(IC,:)-1i*dX(IS,:);

% Handle analysis type
if nargin<=4 || isempty(analysis_type)
    % Default analysis: frequency response
    analysis_type = 'frf';
end
switch lower(analysis_type)
    case {'frf','frequency response'}
        % Frequency response analysis: X = [Q;Om]
        
        % Excitation 'excitation' is the fundamental harmonic of the external
        % forcing
        Fex1 = system.Fex1;
        
        % Excitation frequency
        Om = X(end);
        dOm = zeros(1,length(X)); dOm(end) = 1;
        
        % Scaling of dynamic force equilibrium
        if length(varargin)<2 || isempty(varargin{2})
            fscl = 1;
        else
            fscl = varargin{2};
        end
    otherwise
        error(['Unexpected analysis type ' analysis.type '.']);
end
%% Calculation of the nonlinear forces and the Jacobian
[FnlN,dFnlN] = HB_nonlinear_forces_AFT(QN,dQN,Om,dOm,H,N,...
    system.nonlinear_elements,iN);
%% Assembly of the residual and the Jacobian
smod1 = -Om^2*mmod+1i*Om*dmod+kmod;
dsmod1_dOm = -2*Om*mmod + 1i*dmod;
PHIFex1 = PHI'*Fex1;
QNex = zeros(nN*(H+1),1);
QNex(nN+(1:nN)) = (PHI(iN,:)./smod1)*PHIFex1;
dQNex = zeros(size(QNex,1),size(dOm,2));
dQNex(nN+(1:nN),:) = ...
    (PHI(iN,:).*(-dsmod1_dOm./smod1.^2))*PHIFex1*dOm; % presume dFex = 0
if nN==1
    smod = -Om^2*mmod'*(0:H).^2+1i*Om*dmod'*(0:H)+repmat(kmod',1,H+1);
    dsmod_dOm = -2*Om*mmod'*(0:H).^2 + 1i*dmod'*(0:H);
    HNNinv = diag(reshape(1./((PHI(iN,:).^2)*(1./smod)),nN*(H+1),1));
    dHNN_dOm = diag(reshape((PHI(iN,:).^2)*(-dsmod_dOm./smod.^2),nN*(H+1),1));
else
    % Could probably made more efficient:
    smod = -Om^2*kron((0:H).^2,mmod)+...
        1i*Om*kron((0:H).^1,dmod)+kron((0:H).^0,kmod);
    dsmod_dOm = -2*Om*kron((0:H).^2,mmod) + 1i*kron((0:H).^1,dmod);
    HNN = kron(eye(H+1),PHI(iN,:))*diag(1./smod)*kron(eye(H+1),PHI(iN,:)');
    HNNinv = HNN\eye(nN*(H+1));
    dHNN_dOm = kron(eye(H+1),PHI(iN,:))*diag(-dsmod_dOm./smod.^2)*...
        kron(eye(H+1),PHI(iN,:)');
end
Rc = HNNinv*(QN-QNex) + FnlN;
dRc = HNNinv*(dQN-dQNex) + dFnlN - ...
    HNNinv*( dHNN_dOm * (HNNinv*(QN-QNex)) )*dOm;

% Output Q if requested
if nargout>2
    % Setup excitation vector
    Fex = zeros(n*(H+1),1);
    Fex(n+(1:n)) = Fex1;
    Fnl = zeros(n*(H+1),1);
    iNH = repmat(iN(:),H+1,1) + n*kron((0:H)',ones(nN,1));
    Fnl(iNH) = FnlN;
    smod = -Om^2*kron((0:H).^2,mmod)+...
        1i*Om*kron((0:H).^1,dmod)+kron((0:H).^0,kmod);
    % Could probably made more efficient:
    Q = (kron(eye(H+1),PHI)*diag(1./smod)*kron(eye(H+1),PHI'))*...
        (Fex-Fnl);
end

% Scale dynamic force equilibrium (useful for numerical reasons)
Rc = 1/fscl*(Rc);
dRc = 1/fscl*(dRc);

% Conversion from complex-valued to real-valued residual
R = zeros(size(X,1)-1,1); dR = zeros(size(X,1)-1,size(X,1));
R(I0) = real(Rc(I0)); dR(I0,:) = real(dRc(I0,:));
R(IC) = real(Rc(ID)); dR(IC,:) = real(dRc(ID,:));
R(IS) = -imag(Rc(ID)); dR(IS,:) = -imag(dRc(ID,:));
end
%% Computation of nonlinear forces
function [F,dF] = ...
    HB_nonlinear_forces_AFT(Q,dQ,~,~,H,N,nonlinear_elements,iN)
%% Initialize output

% nonlinear force
F = zeros(size(Q));
dF = zeros(size(F,1),size(dQ,2));

%% Iterate on nonlinear elements
for nl=1:length(nonlinear_elements)
    %% Inverse Discrete Fourier Transform of coordinate associated with
    % this nonlinear element
    
    % Specify time samples along period
    tau = (0:2*pi/N:2*pi-2*pi/N)';
    
    if nonlinear_elements{nl}.islocal
        % Determine force direction associated with nonlinear element
        if size(Q,1)==H+1
            w = 1;
        else
            w = nonlinear_elements{nl}.force_direction(iN);
        end
        W = kron(eye(H+1),w);
        
        % Apply inverse discrete Fourier transform
        H_iDFT = exp(1i*tau*(0:H));
        qnl = real(H_iDFT*(W'*Q));
        dqnl = real(H_iDFT*(W'*dQ));
        %% Evaluate nonlinear force in time domain
        switch lower(nonlinear_elements{nl}.type)
            case 'cubicspring'
                fnl = nonlinear_elements{nl}.stiffness*qnl.^3;
                dfnl = nonlinear_elements{nl}.stiffness*3*...
                    repmat(qnl.^2,1,size(dqnl,2)).*dqnl;
                % For more nonlinear elements, see regular HB_residual.
            otherwise
                error(['Unexpected nonlinear element ' ...
                    nonlinear_elements{nl}.type '.' ...
                    ' For more nonlinear elements, see regular ' ...
                    'HB_residual.']);
        end
        %% Forward Discrete Fourier Transform
        
        % Apply FFT
        Fnlc = fft(fnl(end-N+1:end))/N;
        dFnlc = fft(dfnl(end-N+1:end,:))/N;
        
        % Truncate and convert to half-spectrum notation
        Fnl = [real(Fnlc(1));2*Fnlc(2:H+1)];
        dFnl = [real(dFnlc(1,:));2*dFnlc(2:H+1,:)];
        
        % Store current force into global force vector
        F = F + W*Fnl;
        dF = dF + W*dFnl;
    else % Global nonlinearity
        error('Expecting only local nonlinearities.');
    end
end
end