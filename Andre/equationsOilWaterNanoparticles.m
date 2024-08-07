  function [problem, state] = equationsOilWaterNanoparticles(state0, state, model, ...
                                                      dt, drivingForces, ...
                                                      varargin)
%
% SYNOPSIS:
%   function [problem, state] = equationsOilWaterNanoparticles(state0, state, model, dt, drivingForces, varargin)
%
% DESCRIPTION: 
%   Assemble the linearized equations for an oil-water-nanoparticles
%   system, computing both the residuals and the Jacobians. Returns the result as
%   an instance of the class LinearizedProblem which can be solved using instances
%   of LinearSolverAD.
%
%   A description of the modeling equations can be found in the directory
%   ad-eor/docs.
%
%
% PARAMETERS:
%   state0        - State at previous times-step
%   state         - State at current time-step
%   model         - Model instance
%   dt            - time-step
%   drivingForces - Driving forces (boundary conditions, wells, ...)
%   varargin      - optional parameters
%
% RETURNS:
%   problem - Instance of LinearizedProblem
%   state   - Updated state variable (fluxes, mobilities and more can be
%             stored, the wellSol structure is also updated in case of control switching)
%
% EXAMPLE:
%
% SEE ALSO: LinearizedProblem, LinearSolverAD
%
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

    opt = struct('Verbose'        , mrstVerbose, ... 
                 'reverseMode'    , false      , ... 
                 'velocCompMethod', 'square'   , ...
                 'resOnly'        , false      , ...
                 'iteration'      , -1 );
    opt = merge_options(opt, varargin{:});

    G     = model.G;
    op    = model.operators;
    fluid = model.fluid;
    W     = drivingForces.W;

    % Properties at current timestep
    [p, sW, cs, cs1, cs2] = model.getProps(state, 'pressure', 'water', ...
                                                      'surfactant', ...
                                                      'surfactantdeposition', ...
                                                      'surfactantentrapment');

    % Properties at previous timestep
    [p0, sW0, cs0, cs10, cs20] = model.getProps(state0, 'pressure', 'water', ...
                                                            'surfactant', ...
                                                            'surfactantdeposition', ...
                                                            'surfactantentrapment');

    % Initialize independent variables.
    if ~opt.resOnly
        % ADI variables needed since we are not only computing residuals.
        if ~opt.reverseMode
            [p, sW, cs, cs1, cs2] = initVariablesADI(p, sW, cs, cs1, cs2); 
        else
            % wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
            [p0, sW0, cs0, cs10, cs20] = initVariablesADI(p0, sW0, cs0, cs10, cs20); 
        end
    end

    % We will solve for pressure, water saturation (oil saturation follows via
    % the definition of saturations) and surfactant concentration
    primaryVars = {'pressure', 'sW', 'surfactant', 'surfactantdeposition', 'surfactantentrapment'};

    sO  = 1 - sW;
    sO0 = 1 - sW0;
    sat  = {sW, sO};
    sat0 = {sW0, sO0};
    
    % Update state with AD-variables
    state = model.setProps(state  , {'s', 'pressure', 'surfactant', 'surfactantdeposition', 'surfactantentrapment'}, {sat , p , cs, cs1, cs2});
    state0 = model.setProps(state0, {'s', 'pressure', 'surfactant', 'surfactantdeposition', 'surfactantentrapment'}, {sat0, p0, cs0, cs1, cs2});
    % Set up properties
    state = model.initStateFunctionContainers(state);
    
    % EQUATIONS ---------------------------------------------------------------
    [b, pv]               = model.getProps(state, 'ShrinkageFactors','PoreVolume');
    [b0, pv0]             = model.getProps(state0, 'ShrinkageFactors', 'PoreVolume');
    [phaseFlux, flags]    = model.getProps(state, 'PhaseFlux', 'PhaseUpwindFlag');
    [pressures, mob, rho] = model.getProps(state, 'PhasePressures', 'Mobility', 'Density');
    
    % if isfield(fluid, 'pvMultR')
    %     pvMult =  fluid.pvMultR(p);
    %     pvMult0 = fluid.pvMultR(p0);
    % end
    % pv = pv.*pvMult;
    % pv0 = pv0.*pvMult0;

    [bW, bO]     = deal(b{:});
    [bW0, bO0]   = deal(b0{:});
    [vW, vO]     = deal(phaseFlux{:});
    [upcw, upco] = deal(flags{:});
    [mobW, mobO] = deal(mob{:});
    
    % Upstream weight b factors and multiply by interface fluxes to obtain the
    % fluxes at standard conditions.
    bOvO = op.faceUpstr(upco, bO).*vO;
    bWvW = op.faceUpstr(upcw, bW).*vW;

    % Conservation of mass for water
    water = (1/dt).*( pv.*bW.*sW - pv0.*bW0.*sW0 );
    divWater = op.Div(bWvW);
    water = water + divWater;
    
    % Conservation of mass for oil
    oil = (1/dt).*( pv.*bO.*sO - pv0.*bO0.*sO0 );
    divOil = op.Div(bOvO);
    oil = oil + divOil;
    
    % % Draft ---- Begin:
    % gamma_pt = 1.28;
    % vc = 4.6e-6;
    % gamma_d = 16;
    % gamma_e = 30;
    % % bWvc = op.faceUpstr(upcw, bW).*vc;
    % % csUpstr = op.faceUpstr(upcw,cs);
    % % cs1Upstr = op.faceUpstr(upcw, cs1);
    % 
    % veloc_sq = op.sqVeloc(vW);
    % veloc = veloc_sq.^(1/2);
    % disp(veloc_sq); disp(veloc); disp(vW);
    % 
    % term_cs1_simple = gamma_d.*abs(bW.*veloc).*cs;
    % term_cs2 = gamma_pt.*abs(bW.*veloc).*cs;
    % term_cs1_compl = gamma_e.*abs(bW.*(veloc-vc)).*cs1;
    % 
    % if vW <= vc
    %     RHS_cs1 = term_cs1_simple;
    %     RHS = term_cs1_simple + term_cs2;
    % else
    %     RHS_cs1 = term_cs1_simple - term_cs1_compl;
    %     RHS = term_cs1_simple - term_cs1_compl + term_cs2;
    % end
    % 
    % 
    % vSft   = op.faceUpstr(upcw, cs).*vW;
    % bWvSft = op.faceUpstr(upcw, bW).*vSft;
    % D = 5.6e-8;
    % pvbWsWD= op.faceUpstr(upcw, pv.*bW.*sW.*D);
    % pvbWsWDcs = pvbWsWD.*op.Grad(cs);
    % 
    % surfactant = (1/dt).*(pv.*bW.*sW.*cs - pv0.*bW0.*sW0.*cs0);
    % divSurfactant = op.Div(bWvSft-pvbWsWDcs);
    % surfactant = surfactant + divSurfactant;
    % 
    % % 2 equations with 2 unknowns i.e., cs and cs1:
    % eq_1 = surfactant - RHS;
    % eq_2 = (1/dt).*(cs1 -cs10) - RHS_cs1;
    % 
    % eq_3 = (1/dt).*(cs2-cs20) - term_cs2;
    % 
    % % Draft ---- END

    % EOR EQUATIONS ---------------------------------------------------------------
    % Set parameters
    vc = 4.6e-6; % critical velocity for the water phase
    gamma_d = 16; % rate coefficients for surface retention of the nanoparticles in the water phase
    gamma_e = 30; % rate coefficients for entrainment of the nanoparticles in the water phase
    gamma_pt = 1.28; % pore throat blocking constant
    veloc_sq = op.sqVeloc(vW);
    veloc = veloc_sq.^(1/2);
    %    We use a velocity reconstruction of the type  v_c  = 1/V * sum_{f} (x_f -x_c)u_f
    % where : v_c : Approximated alue of the velocity at the cell center
    % V   : Cell volume, f   : Face, x_f : Centroid of the face f
    % x_c : Centroid of the cell, u_f : flux at the face f

    % RHS of equations:
    if vW <= vc
        derDeposition = gamma_d.*abs(bW.*veloc).*cs;
        derEntrapment = gamma_pt.*abs(bW.*veloc).*cs;
    else
        derDeposition = (gamma_d.*abs(bW.*veloc).*cs) + (gamma_e.*abs(bW.*(veloc-vc)).*cs1);
        derEntrapment = gamma_pt.*abs(bW.*veloc).*cs;
    end
    R = derDeposition + derEntrapment; 

    % Deposition equation:
    deposition = (1/dt).*(cs1 - cs10) - derDeposition;

    % entrapment equation
    entrapment = (1/dt).*(cs2 - cs20) - derEntrapment;

    % Conservation of surfactant in water:
    vSft   = op.faceUpstr(upcw, cs).*vW;
    bWvSft = op.faceUpstr(upcw, bW).*vSft;
    D = 5.6e-8; % *ones(G.cells.num, 1);
    % D = initVariablesADI(D);
    pvbWsWD= op.faceUpstr(upcw, pv.*bW.*sW.*D);
    pvbWsWDcs = pvbWsWD.*op.Grad(cs);

    surfactant    = (1/dt).*(pv.*bW.*sW.*cs - pv0.*bW0.*sW0.*cs0);
    divSurfactant = op.Div(bWvSft-pvbWsWDcs); %Not sure about porosity and faceUpstr
    surfactant = surfactant + divSurfactant - R;
    %-------------------------------------------------------------------------------------------
    % Computation of adsoprtion term
    %poro = model.rock.poro
    %ads  = model.getProp(state , 'SurfactantAdsorption');
    %ads0 = model.getProp(state0, 'SurfactantAdsorption');
    %ads_term = fluid.rhoRSft.*((1-poro)./poro).*(ads - ads0);
    %
    % surfactant    = (1/dt).*(pv.*bW.*sW.*cs - pv0.*bW0.*sW0.*cs0); % + (op.pv/dt).*ads_term;
    % divSurfactant = op.Div(bWvSft);
    % surfactant = surfactant + divSurfactant;
    % fprintf('after: %f\n %f\n %f\n %f\n %f\n %f\n %f\n %f\n %f\n %f\n \n', size(D), size(bW),size(sW),size(pv), size(cs),size(bW0),size(sW0),size(pv0), size(cs0), size(dt))

    if ~opt.resOnly
        epsilon = 1.e-8;
        % the first way is based on the diagonal values of the resulting
        % Jacobian matrix
        eps = sqrt(epsilon)*mean(abs(diag(surfactant.jac{3})));
        % sometimes there is no water in the whole domain
        if (eps == 0.)
            eps = epsilon;
        end
        % bad marks the cells prolematic in evaluating Jacobian
        bad = abs(diag(surfactant.jac{3})) < eps;
        % the other way is to choose based on the water saturation
        surfactant(bad) = cs(bad);
    end

    % % Draft ---- Begin:
    % 
    % eqs      = {water   , oil   , eq_1, eq_2, eq_3};
    % names    = {'water' , 'oil' , 'surfactant', 'surfactantdeposition', 'surfactantentrapment'};
    % types    = {'cell'  , 'cell', 'cell', 'cell', 'cell'};
    % components = {cs, cs1, cs2};
    % 
    % % Draft ---- END

    eqs      = {water   , oil   , surfactant, deposition, entrapment};
    names    = {'water' , 'oil' , 'surfactant', 'surfactantdeposition', 'surfactantentrapment'};
    types    = {'cell'  , 'cell', 'cell', 'cell', 'cell'};
    components = {cs, cs1, cs2};
     
    [eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, ...
                                                   state, pressures, sat, mob, ...
                                                   rho, {}, components, ...
                                                   drivingForces);
    % % Finally, add in and setup well equations
    % [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, ...
    %                                                   types, wellSol0, wellSol, ...
    %                                                   wellVars, wellMap, p, ...
    %                                                   mob, rho, {}, components, ...
    %                                                   dt, opt);

    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end


