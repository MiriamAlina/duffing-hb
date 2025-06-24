%========================================================================
% DESCRIPTION:
% Class representing a reduced order model derived from an object of class
% 'MechanicalSystem' defined in NLvib. Four types of component mode 
% synthesis methods are implemented, specified by the 'type' input:
%       'CB'    (conventional) Craig-Bampton method,
%       'CBMB'  massless-boundary Craig-Bampton method,
%       'RU'    Rubin method, and
%       'MCN'   MacNeal method.
% All those methods are defined in [2]. In each of those methods, a set of
% coordinates of the parent finite element model is to be specified as 
% contact boundary coordinates via the index vector 'IB', and a
% number of 'nm' normal modes is retained (e.g. for fixed-boundary in the 
% 'CB' case, and for free-boundary in the 'RU' case). The contact boundary
% coordinates are retained in the reduced model thanks to the use of 
% appropriate static modes (e.g., static constraint modes in the 'CB' case, 
% residual flexibility attachment modes in the 'RU' case). Static and 
% (dynamic) normal modes form the component modes, which are stored as
% columns of the matrix 'T'. The reduced stiffness, damping and mass
% matrices are 'M', 'D', and 'K', respectively.
% 
% The primary purpose of this class is model reduction in the context of
% dynamic contact simulation with the tool NLstep [1]. This is where the
% massless-boundary methods 'CBMB' and 'MCN' are particularly useful [2].
% Note that 'CB' and 'RU' are Galerkin-consistent methods in the sense 
% that the reduced matrices are obtained by the inner product with the 
% matrix 'T'. This does not hold for the mass-deficient methods 'CBMB' and
% 'MCN'.
%
% [1] https://github.com/maltekrack/
% [2] https://doi.org/10.1016/j.compstruc.2021.106698
% 
% Mandatory variables:
%   model               MechanicalSystem object
%   IB                  index vector of boundary coordinates
%   nm                  integer specifying number or retained normal modes
%   type                string for the model order reduction method
%                       available options: 'CB','CBMB','RU','MCN'
%========================================================================
% This file is part of NLvib.
%
% If you use NLvib, please refer to the book:
%   M. Krack, J. Gross: Harmonic Balance for Nonlinear Vibration
%   Problems. Springer, 2019. https://doi.org/10.1007/978-3-030-14023-6.
%
% COPYRIGHT AND LICENSING:
% NLvib Copyright (C) 2025
%   Malte Krack (malte.krack@ila.uni-stuttgart.de)
%   Johann Gross (johann.gross@ila.uni-stuttgart.de)
%   University of Stuttgart
% This program comes with ABSOLUTELY NO WARRANTY.
% NLvib is free software, you can redistribute and/or modify it under the
% GNU General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% For details on license and warranty, see http://www.gnu.org/licenses
% or gpl-3.0.txt.
%========================================================================

classdef CMS_ROM < handle

    properties
        type                % type (string) of model reduction method
        n                   % total number of generalized coordinates
        nb                  % number of contact boundary coordinates
        T                   % matrix of component modes
        K                   % reduced stiffness matrix
        M                   % reduced mass matrix
        D                   % reduced damping matrix
        Fex1                % excitation force vector
        Tresp               % response vector (specific row of 'T')
        nonlinear_elements  % structure defining nonlinear elements
        rbmTol=.1;          % threshold for modal angular frequency to be 
                            % considered as rigid body mode
    end

    properties (Constant)
        MBtypes = ...       % massless-boundary model reduction types
            {'CBMB','MCN'};
    end

    methods
        function obj = CMS_ROM(model,IB,nm,type,rbmTol)
            %% Constructor
            obj = obj@handle;

            % Check input
            if nargin==0
                return
            end
            if ~isa(model,'MechanicalSystem')
                error(['First input must be an object of ' ...
                    'MechanicalSystem class.']);
            end

            % If no threshold has been specified for recognizing rigid body
            % modes, use default value
            if nargin<5 || isempty(rbmTol)
                rbmTol = obj.rbmTol;
            end

            % Store type and build reduced order model
            obj.type = type;
            [Mr,Kr,Tr] = CMS_ROM.buildROM(model.M,model.K,...
                IB,nm,type,rbmTol);

            % Store reduction and reduced mass and stiffness matrices
            obj.T = Tr;
            obj.M = Mr;
            obj.K = Kr;

            % Project damping matrix and excitation force vector
            obj.D = Tr'*model.D*Tr;
            obj.Fex1 = Tr'*model.Fex1;

            % Determine total number of generalized coordinates and number
            % of contact boundary coordinates (=number of static modes)
            obj.n = length(Mr);
            obj.nb = length(IB);

            % Adopt nonlinear elements from parent model, check and project
            obj.nonlinear_elements = model.nonlinear_elements;
            check_nonlinearities(obj);
            transform_nonlinear_force_directions(obj,Tr);
        end
        function apply_transform(obj,Tr)
            obj.T = obj.T*Tr;
            obj.M = Tr'*obj.M*Tr;
            obj.D = Tr'*obj.D*Tr;
            obj.K = Tr'*obj.K*Tr;
            obj.Fex1 = Tr'*obj.Fex1;
            obj.Tresp = obj.Tresp*Tr;
            transform_nonlinear_force_directions(obj,Tr);
        end
        function check_nonlinearities(obj)
            for nl = 1:length(obj.nonlinear_elements)

                % check 3D-contact-specific properties
                if strncmp(obj.nonlinear_elements{nl}.type,'3D',2)

                    % Contact is a local nonlinearity, and friction makes
                    % it hysteretic in the sense of the tool
                    obj.nonlinear_elements{nl}.islocal = 1;
                    obj.nonlinear_elements{nl}.ishysteretic = 1;
                    
                    % Exactly one specificiation is needed, either
                    % 'preload' or 'imposedGap'.
                    if isfield(obj.nonlinear_elements{nl},'preload') ...
                            && isfield(obj.nonlinear_elements{nl},'imposedGap')
                        error('Cannot have preload and imposed gap simulstaneously.');
                    elseif ~isfield(obj.nonlinear_elements{nl},'preload') ...
                            && ~isfield(obj.nonlinear_elements{nl},'imposedGap')
                        error('Either preload or imposed gap must be specified.');
                    end
                elseif strncmp(obj.nonlinear_elements{nl}.type,'1D',2) ||...
                        strncmp(obj.nonlinear_elements{nl}.type,'friction',8) ||...
                        strncmp(obj.nonlinear_elements{nl}.type,'unilateral',10)

                    % Contact is a local nonlinearity, and friction makes
                    % it hysteretic in the sense of the tool
                    obj.nonlinear_elements{nl}.islocal = 1;
                    
                    % Exactly one specificiation is needed, either
                    % 'friction_limit_force' or 'imposedGap'.
                    if isfield(obj.nonlinear_elements{nl},'friction_limit_force') ...
                            && isfield(obj.nonlinear_elements{nl},'imposedGap')
                        error('Cannot have preload and imposed gap simulstaneously.');
                    elseif ~isfield(obj.nonlinear_elements{nl},'friction_limit_force') ...
                            && ~isfield(obj.nonlinear_elements{nl},'imposedGap')
                        error('Either preload or imposed gap must be specified.');
                    end

                    if isfield(obj.nonlinear_elements{nl},'friction_limit_force')
                        obj.nonlinear_elements{nl}.ishysteretic = 1;
                    else
                        obj.nonlinear_elements{nl}.ishysteretic = 0;
                    end
                end
            end
        end
        function transform_nonlinear_force_directions(obj,Tr)
            % Project force direction of nonlinear elements, if specified
            for nl = 1:length(obj.nonlinear_elements)
                if isfield(obj.nonlinear_elements{nl},'force_direction')
                    obj.nonlinear_elements{nl}.force_direction = ...
                        Tr'*obj.nonlinear_elements{nl}.force_direction;
                end
            end
        end
        function isMB = isMasslessBoundaryROM(obj)
            isMB = any(strcmpi(obj.MBtypes,obj.type));
        end
    end
    methods (Static)
        function [Mr,Kr,Tr] = buildROM(M,K,IB,nm,type,rbmTol)
            % Check input
            n = length(M);
            if size(M,1)~=n || size(M,2)~=n || ...
                    size(K,1)~=n || size(K,2)~=n
                error('M and K must be n x n matrices.');
            end
            % If contact boundary coordinates are specified as logical, 
            % convert to index vector
            if islogical(IB)
                IB = find(IB);
            end
            nb = length(IB);
            if numel(nm)~=1 || ~nm>0 || ~mod(nm,1)==0
                error('Number of dynamic modes must be positive integer.');
            end
            if (nb+nm)>n
                error(['Number of contact boundary coordinates ' ...
                    'plus number of dynamic modes cannot be larger ' ...
                    'than the number of coordinates of the parent model.']);
            end
            if numel(rbmTol)~=1 || ~rbmTol>0
                error(['Rigid body frequency threshold must be ' ...
                    'positive scalar.']);
            end

            % Inform user about reduction parameters
            disp(['**INFO: Parent model has ' num2str(n) ...
                ' degrees of freedom. Performing ' type ...
                '-reduction with ' num2str(length(IB)) ' static and ' ...
                num2str(nm) ' dynamic modes. ']);

            % First we rearrange so that contact boundary coordinates are
            % at the start, and the remaining coordinates follow.
            II = setdiff(1:n,IB)';
            order = [IB;II];
            Tsort = sparse(1:n,order,1)';
            M = Tsort'*M*Tsort;
            K = Tsort'*K*Tsort;

            % Thanks to the above rearragement, we now have trivial index
            % vectors
            IB = 1:nb;
            II = (nb+1):n;
            switch upper(type)
                case {'CB','CBMB'}
                    % Compute static constraint modes
                    PSI = -K(II,II)\K(II,IB);

                    % Compute 'nm' lowest-frequency fixed-interface 
                    % normal modes
                    [THETA,OM2] = eigs(K(II,II),M(II,II),nm,'sm');
                    om = sqrt(diag(OM2));
                    % Sort and normalize w.r.t. mass matrix
                    [om,ind] = sort(om); THETA = THETA(:,ind);
                    THETA = THETA./repmat(sqrt(diag(THETA'*M(II,II)*THETA))',...
                        size(THETA,1),1);

                    % Display lowest/highest retained frequency
                    fprintf(['**INFO: lowest frequency among retained ' ...
                        'normal modes: EFmax=%.2f [Hz]\n'], min(om)/2/pi);
                    fprintf(['**INFO: highest frequency among retained ' ...
                        'normal modes: EFmax=%.2f [Hz]\n'], max(om)/2/pi);
                    
                    % Set up reduced mass and stiffness matrices
                    switch upper(type)
                        case 'CB'
                            Mbb = M(IB,IB) + M(IB,II)*PSI + ...
                                PSI'*M(II,IB) + PSI'*M(II,II)*PSI;
                            Mbi = (M(IB,II) + PSI'*M(II,II))*THETA;
                            Kbb = K(IB,IB) + K(IB,II)*PSI + ...
                                PSI'*K(II,IB) + PSI'*K(II,II)*PSI;
                            Mr = [Mbb Mbi; Mbi' eye(nm)];
                            Kr = blkdiag(Kbb,diag(om.^2));
                        case 'CBMB'
                            alpha = THETA'*(M(II,IB) + M(II,II)*PSI);
                            PSI = PSI-THETA*alpha;
                            Kbb = K(IB,IB) + K(IB,II)*PSI + ...
                                PSI'*K(II,IB) + PSI'*K(II,II)*PSI;
                            Kbi = (K(IB,II) + PSI'*K(II,II))*THETA;
                            Mr = blkdiag(zeros(nb),eye(nm));
                            Kr = [Kbb Kbi; Kbi' diag(om.^2)];
                    end
                    
                    % Set up matrix of component modes
                    Tr = zeros(n,nb+nm);
                    Tr(IB,1:nb) = eye(nb);
                    Tr(II,1:nb) = PSI;
                    Tr(II,nb+(1:nm)) = THETA;
                case {'MCN','RU'}
                    % Compute 'nm' lowest-frequency free-interface 
                    % normal modes
                    [PHI,OM2] = eigs(K,M,nm,'sm');
                    om = sqrt(diag(OM2));
                    % Sort and normalize w.r.t. mass matrix
                    [~,ind] = sort(om); PHI = PHI(:,ind);
                    PHI = PHI./repmat(sqrt(diag(PHI'*M*PHI))',size(PHI,1),1);

                    % Display lowest/highest retained frequency
                    fprintf(['**INFO: lowest frequency among retained ' ...
                        'normal modes: EFmax=%.2f [Hz]\n'], min(om)/2/pi);
                    fprintf(['**INFO: highest frequency among retained ' ...
                        'normal modes: EFmax=%.2f [Hz]\n'], max(om)/2/pi);

                    % Account for potential rigid body modes via inertia
                    % relief
                    isRBM = (om<rbmTol);    
                    if any(isRBM)
                        P = eye(n)-M*PHI(:,isRBM)*PHI(:,isRBM)';
                    else
                        P = speye(n);
                    end
                    % Calculate flexibility sub-matrix associated with
                    % contact boundary
                    Lb = P*speye(n,nb);
                    F = K\Lb;
                    % Calculate corresponding residual flexibility matrix
                    dF = F - PHI(:,~isRBM)*diag(1./om(~isRBM).^2)*...
                        PHI(IB,~isRBM)';

                    % Set up matrix of component modes
                    T1 = zeros(n,nb+nm);
                    T1(:,1:nb) = P'*dF;
                    T1(:,nb+(1:nm)) = PHI;
                    dFbb = Lb'*dF;
                    dKbb = dFbb\eye(length(dFbb));
                    T2 = zeros(size(T1,2),nb+nm);
                    T2(:,1:nb) = [dKbb;zeros(nm,nb)];
                    T2(:,nb+(1:nm)) = [-dKbb*PHI(IB,:);eye(nm)];
                    Tr = T1*T2;

                    % Set up reduced mass and stiffness matrix
                    switch upper(type)
                        case 'MCN'
                            Mr = blkdiag(zeros(nb), eye(nm));
                            Kr = [dKbb -dKbb*PHI(IB,:); ...
                                -PHI(IB,:)'*dKbb ...
                                diag(om.^2) + PHI(IB,:)'*dKbb*PHI(IB,:)];
                        case 'RU'
                            Mr = Tr'*M*Tr;
                            Kr = Tr'*K*Tr;
                    end
                otherwise
                    error(['Unknown reduction type specifier ' type '.']);
            end

            % The matrix of component modes lives in the initial space of
            % the parent model. Hence we have to apply the sorting, too.
            Tr = Tsort * Tr;
        end
    end
end
