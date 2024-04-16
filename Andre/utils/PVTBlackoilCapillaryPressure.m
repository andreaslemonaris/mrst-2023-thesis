classdef PVTBlackoilCapillaryPressure < StateFunction
    % Dummy PVT capillary pressure
    properties
    end
    
    methods
        function gp = PVTBlackoilCapillaryPressure(model, varargin)
            gp@StateFunction(model, varargin{:});
            assert(isfield(model.fluid, 'pcOW'), ...
                'No capillary pressure in the fluid object');
            gp = gp.dependsOn('CapillaryPressure', 'FlowPropertyFunctions');
            gp.label = '\Pc_OW';
        end
        function pv = evaluateOnDomain(prop, model, state) %#ok
            % Dummy PVT capillary pressure
            p_c = prop.getEvaluatedExternals(model, state, 'BlackOilCapillaryPressure');
            p_c_PVT = p_c;
        end
    end
end