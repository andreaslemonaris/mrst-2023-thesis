classdef PVTBaseRelativePermeability < StateFunction
    % Dummy PVT relative permeability
    properties
    end
    
    methods
        function gp = PVTBaseRelativePermeability(model, varargin)
            gp@StateFunction(model, varargin{:});
            assert(isfield(model.fluid, 'krW'), ...
                'Absence of krW in the fluid object');
             assert(isfield(model.fluid, 'krO'), ...
                'Absence of krO in the fluid object');
            gp = gp.dependsOn('RelativePermeability', 'FlowPropertyFunctions');
            gp.label = '\kr_OW';
        end
        function pv = evaluateOnDomain(prop, model, state) %#ok
            % Dummy PVT capillary pressure
            k_r = prop.getEvaluatedExternals(model, state, 'BaseRelativePermeability');
            k_r_PVT = k_r;
        end
    end
end