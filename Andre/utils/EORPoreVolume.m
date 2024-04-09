classdef EORPoreVolume < PoreVolume
    % Effective pore-volume. Accounts for effect the reduction of porosity
    % that nanoparticles deposition on the pore surfaces and/or blocking of pore throats induces.
    properties
    end
    
    methods
        function gp = EORPoreVolume(model, varargin)
            gp@PoreVolume(model, varargin{:});
            gp = gp.dependsOn({'pressure'}, 'state');
            gp = gp.dependsOn({'surfactantdeposition'}, 'state');
            gp = gp.dependsOn({'surfactantentrapment'}, 'state');
            assert(isfield(model.water && model.surfactant), 'surfactant is missing'); %check mechanism
 %           assert(isfield(model.fluid, 'pvMultR'), 'pvMultR missing from fluid.');
        end
        function pv = evaluateOnDomain(prop, model, state)
            % Get effective pore-volume, accounting for rock-compressibility
            pv = evaluateOnDomain@PoreVolume(prop, model, state);
            c1 = model.getProps(state, 'surfactantdeposition');
            c2 = model.getProps(state, 'surfactantentrapment');

 %           p = model.getProp(state, 'pressure');
 %           pvMult = prop.evaluateFluid(model, 'pvMultR', p);
 %           pv = pv.*pvMult;
            pv = pv - (c1 + c2);
                
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
