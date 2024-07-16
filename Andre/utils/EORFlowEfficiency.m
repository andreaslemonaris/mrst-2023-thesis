classdef EORFlowEfficiency < StateFunction
    % Flow efficiency factor. Expressing the fraction of unplugged pores available for flow.
    properties
        gamma_f % Coefficient of flow efficiency 
    end
    
    methods
        function gp = EORFlowEfficiency(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'surfactantentrapment'}, 'state');
            % assert(isfield(model.surfactant), 'surfactant is missing'); %check mechanism
            gp.gamma_f = 0.01;
        end
        function f = evaluateOnDomain(prop, model, state) %#ok
            % Get flow efficiency
            c2 = model.getProps(state, 'nanoparticlesentrapment');
            gammaf = prop.gamma_f;

            f = 1 - gammaf.*c2;      
        end
    end
end

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}