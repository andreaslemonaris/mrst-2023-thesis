%% Single-phase slightly compressible fluid (water)
%% Input concentration for formation and injected waters is in mole/kilogram
%% Cartesian grid model
clc, clear, close all
tic
%% Define and descretize the domain
[nx,ny,nz] = deal(50, 1, 1);
[Lx,Ly,Lz] = deal(500*ft, 1*ft, 1*ft);
G = cartGrid([nx, ny, nz], [Lx, Ly, Lz]);
G = computeGeometry(G);

plotGrid(G); view(3); axis tight

%% Define rock model
rock   = makeRock(G, 30*milli*darcy, 0.3);
cr     = 5e-6/psia;
Pref   = 3000*psia;
Dref   = 0.5*ft;
PVref  = poreVolume(G, rock);
PVfunc = @(P) PVref.*exp(cr*(P - Pref));

P = linspace(1500*psia,3000*psia,50);
plot(P/psia, PVref(1)*exp(cr*(P-Pref))/ft^3,'LineWidth',2);
xlabel('Pressure [psia]')
ylabel('Single Gridblock Pore Volume [ft^3]')

%% Define model for an slightly compressible fluid
mu     = 1*centi*poise;
cFluid = 3e-6/psia;
rhoRef = 62*pound/ft^3;    % reservoir fluid density at the reference conditions
rhoFunc = @(P) rhoRef*exp(cFluid*(P - Pref));

plot(P/psia,rhoFunc(P)/(pound/ft^3),'LineWidth',2);
xlabel('Pressure [psia]')
ylabel('Fluid Density [lbm/ft^3]')

%% Specify boundary conditions
nc = G.cells.num;
W = addWell([], G, rock,  1, 'Name', 'Injection Well' , 'Radius', 0.5*ft);
W = addWell(W , G, rock, nc, 'Name', 'Production Well', 'Radius', 0.5*ft);
injCell  = W(1).cells;
prodCell = W(2).cells;
PIprod = W(2).WI;    % well productivity indices
dzProd = W(2).dZ;    % perforated cell depth relative to bottom-hole
% Internal boundary conditions
qInj = 0.1*stb/day;
BHPprod = Pref;
primNames = {'H', 'CO3', 'Ca'};
influxConc = [1, 1, 1]*mol/kilogram;

gravity reset on, g = norm(gravity);
Pwellbore  = BHPprod + rhoFunc(BHPprod)*g.*dzProd; % wellbore pressures at perforated points

%% Set initial pressure and concentration
PG = rhoFunc(Pref)*g;
Pinit = Pref + PG.*(G.cells.centroids(:,3) - Dref);

Cinit = [0.5, 0.5, 0.5]*mol/kilogram;
Cinit = repmat(Cinit, nc, 1);

%% Compute transmissibilities
N  = double(G.faces.neighbors);             % Map: face -> cell
intInx = all(N ~= 0, 2);                    % Interior faces
N  = N(intInx, :);                          % Interior neighbors
hT = computeTrans(G, rock);                 % Half-transmissibilities
cf = G.cells.faces(:,1);                    % Map: cell -> face
nf = G.faces.num;                           % Total number of faces
T  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]); % Harmonic average
T  = T(intInx);                             % Restricted to interior

%% Define discrete operators to desretize the mathematical equations
n = size(N,1);
C = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[-1 1], n, G.cells.num);
grad    = @(x)C*x;
div     = @(x)-C'*x;
avg     = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));
upwind  = @(x,flag) x(N(:,1)).*double(flag)+x(N(:,2)).*double(~flag);

%% Initialize for solution loop
[P_AD, cH_AD, cCO3_AD, cCa_AD] = initVariablesADI(Pinit, Cinit(:,1), Cinit(:,2), Cinit(:,3));
C_AD = cell(1, numel(primNames));
[C_AD{1}, C_AD{2}, C_AD{3}] = deal(cH_AD, cCO3_AD, cCa_AD);
[pIx, cHIx, cCO3Ix, cCaIx] = deal(1:nc, (nc+1):(2*nc), (2*nc+1):(3*nc), (3*nc+1):(4*nc));

numSteps = 50;                  % number of time-steps
totTime  = 150*day;             % total simulation time
dt       = totTime / numSteps;  % constant time step
tol      = 1e-5;                % Newton tolerance
maxIter  = 10;                  % max number of Newton-Raphson iterations

sol = repmat(struct('time',[],'pressure',[], 'C', []),[numSteps+1,1]);
sol(1)  = struct('time', 0, 'pressure', value(P_AD)/psia, 'C', [value(cH_AD), value(cCO3_AD), value(cCa_AD)]);

%% Main loop
t = 0; step = 0;
hwb = waitbar(t,'Simulation ...');
while t < totTime
    t = t + dt;
    step = step + 1;
    fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
        step, convertTo(t - dt, day), convertTo(t, day));

    % Newton loop
    resNorm = 1e+99;
    P0  = value(P_AD); % Previous step pressure
    cH0  = value(C_AD{1});
    cCO30  = value(C_AD{2});
    cCa0  = value(C_AD{3});
    C0 = [cH0, cCO30, cCa0];
    iterCounter = 0;
    while (resNorm > tol) && (iterCounter <= maxIter)

        %% Define pressure equation
        [rho_AD, rho0] = deal(rhoFunc(P_AD), rhoFunc(P0));
        [PV_AD , PV0]  = deal(PVfunc(P_AD), PVfunc(P0));
        gradz  = grad(G.cells.centroids(:,3));
        v = -(T/mu).*(grad(P_AD) - avg(rho_AD)*g.*gradz);
        pressEq = (1/dt)*(PV_AD.*rho_AD - PV0.*rho0) + div(avg(rho_AD).*v);

        % Adding boundary conditions
        qProd  = -PIprod.*(1/mu).*(P_AD(prodCell) - Pwellbore);
        % The BC of the production well can be set in 2 different ways, either using the BHP directly or using the equivalent production flowrate.
        pressEq(prodCell) = pressEq(prodCell) - qProd*rho_AD(prodCell); % Method 1
        %       pressEq(prodCell) = P_AD(prodCell) - BHPprod;      % Method 2

        pressEq(injCell)  = pressEq(injCell) - qInj*rhoRef;

        %% Define the mass transport equation of the desired species
        % water upstream-index
        upcw = v > 0;
        transEq = cell(1,numel(primNames));
        for i = 1:numel(primNames)
            accumTerm  = (1/dt).*(PV_AD.*rho_AD.*C_AD{i} - PV0.*rho0.*C0(:,i));
            % Flux term without dispersion
%             fluxTerm   = div(upwind(C_AD{i}.*rho_AD, upcw).*v);
            % Flux term with dispersion
            D = [10^(-7), 10^(-8), 10^(-9)]*meter^2/second*ones;
            D = D(i);
            fluxTerm   = div(upwind(C_AD{i}.*rho_AD, upcw).*v - avg(rock.poro.*rho_AD).*D.*grad(C_AD{i}));

            transEq{i} = accumTerm + fluxTerm;

            % Add boundary conditions
%             transEq{i}(injCell)  = C_AD{i}(injCell) - influxConc(i);   % Method 1
            transEq{i}(injCell) = transEq{i}(injCell) - qInj*rhoRef*influxConc(i);   % Method 2
            transEq{i}(prodCell) = transEq{i}(prodCell) - qProd.*rho_AD(prodCell)*C_AD{i}(prodCell);
        end

        %% Collect and concatenate all equations (i.e., assemble and linearize system)
        eqs = {pressEq, transEq{1}, transEq{2}, transEq{3}};
        eq  = cat(eqs{:});
        J   = eq.jac{1};  % Jacobian
        res = eq.val;     % residual
        increment = -(J \ res); % Newton update


        %% Update variables
        P_AD.val   = P_AD.val   + increment(pIx);
        C_AD{1}.val   = C_AD{1}.val   + increment(cHIx);
        C_AD{2}.val   = C_AD{2}.val   + increment(cCO3Ix);
        C_AD{3}.val   = C_AD{3}.val   + increment(cCaIx);

        resNorm = norm(res);
        iterCounter   = iterCounter + 1;
        fprintf('  Iteration %3d:  Res = %.4e\n', iterCounter, resNorm);
    end

    if iterCounter > maxIter
        error('Newton solver did not converge')
    else % Store solution
        sol(step+1)  = struct('time', t, 'pressure', value(P_AD)/psia, 'C', [value(C_AD{1}), value(C_AD{2}), value(C_AD{3})]);
        waitbar(t/totTime,hwb);
    end
end
runTime = toc;
fprintf('\nSimulation is completed\n')
close(hwb);

%% For a 1-dimensional domain
clf
for i = 1:numel(sol)
    for j = 1:numel(primNames)
        subplot(1, numel(primNames), j)

        plot(1:nc, sol(i).C(:,j), 'LineWidth', 1)
        title(['Tracer Concentration Profile (', 'Time Step: ', num2str(i),')'])
        xlim([1,nc])
        ylim([Cinit(1,j), influxConc(j)])
        xlabel('Gridblock Index')
        ylabel('Tracer Concentration [mole/kilogram]')

    end
    pause(0.25)
end

%% mrst-gui
mrstModule add mrst-gui
plotToolbar(G, sol)
view(3)
colorbar
colormap(jet)
shading faceted
axis tight off
