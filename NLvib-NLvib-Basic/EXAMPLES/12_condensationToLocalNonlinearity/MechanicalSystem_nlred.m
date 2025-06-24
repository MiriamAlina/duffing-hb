%========================================================================
% DESCRIPTION: 
% This class extends the MechanicalSystem class by a few additional
% properties which are required for spectral decomposition of the dynamic
% compliance matrix.
%========================================================================
classdef MechanicalSystem_nlred < MechanicalSystem
    properties
        iN
        mmod
        kmod
        dmod
        PHI
    end
    
    methods
        function obj = MechanicalSystem_nlred(M,D,K,...
                nonlinear_elements,Fex1)
            %% Constructor
            obj = obj@MechanicalSystem(M,D,K,nonlinear_elements,Fex1);
        end
    end
end