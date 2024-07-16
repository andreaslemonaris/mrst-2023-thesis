classdef EORAbsolutePermReduction < StateFunction
% Absolute permeability reduction. Accounts for effect the reduction of
% absolute permeability that nanoparticles deposition on the pore surfaces
% and/or blocking of pore throats induces.
    
    properties
        % constant for fluid seepage allowed by the plugged pores
        k_f
        %exponent l of the equation -typical value of the range from 2.5 to
        %3.5
        exponent_l 
    end

    methods
        function gp = EORAbsolutePermReduction(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PoreVolume', 'FlowEfficiencyFactor', 'BasePoreVolume'});
            gp = gp.dependsOn({'surfactantdeposition'}, 'state');
            gp = gp.dependsOn({'surfactantentrapment'}, 'state');
            % assert(isfield(model.water && model.surfactant), 'surfactant is missing'); %check mechanism
            gp = gp.dependsOn('Transmissibility', 'FlowDiscretization');
            gp.k_f = 0.6;
            gp.exponent_l = 3;
            gp.label = '\perm_\alpha';
        end

        function permRed = evaluateOnDomain(prop, model, state)
            [pv, f] = prop.getEvaluatedDependencies(state, 'PoreVolume', 'FlowEfficiencyFactor');
            T = prop.getEvaluatedExternals(model, state, 'Transmissibility');
            kf = prop.k_f;
            l = prop.exponent_l;
            poro = prop.getEvaluatedDependencies(state, 'BasePoreVolume');
            s = model.operators;
            % poro = model.rock.poro;
            red = ((1-f).*kf + f.*(pv./poro)).^l;
            red = s.faceAvg(red);
            permRed = T.*red;
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
