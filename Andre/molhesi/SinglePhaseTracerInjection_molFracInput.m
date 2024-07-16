%% Single-phase slightly compressible fluid (water)
%% Input concentration for formation and injected waters is in mole fraction
%% Cartesian grid model
clc, clear, close all

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
qInj    = 0.1*stb/day;
BHPprod = Pref;
MW_H2O = 18.01528*gram/mol;
MW_Ca  = 40.078*gram/mol;   % If solute is Calcium
xInflux = 0.05;   % mole fraction of solute in injected water
Xinflux = xInflux*MW_Ca/(1*MW_H2O + 0*MW_Ca);   % mass fraction of slute (for the denominator, the amount of solute is ignored)

gravity reset on, g = norm(gravity);
Pperf  = BHPprod + rhoFunc(BHPprod)*g.*dzProd; % wellbore pressures at perforated points

%% Set initial pressure and concentration
PG = rhoFunc(Pref)*g;
Pinit = Pref + PG.*(G.cells.centroids(:,3) - Dref);

xInit = 0.01*ones(nc,1);   % mole fraction of solute in formation water

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
[P_AD, x_AD] = initVariablesADI(Pinit, xInit);
[pIx, xIx] = deal(1:nc, (nc+1):(2*nc));

numSteps = 50;                  % number of time-steps
totTime  = 150*day;             % total simulation time
dt       = totTime / numSteps;  % constant time step
tol      = 1e-5;                % Newton tolerance
maxIter  = 10;                  % max number of Newton-Raphson iterations

sol = repmat(struct('time',[],'pressure',[], 'x', []),[numSteps+1,1]);
sol(1)  = struct('time', 0, 'pressure', value(P_AD)/psia, 'x', value(x_AD));

%% Main loop
t = 0; step = 0;
hwb = waitbar(t,'Simulation ...');
while t < totTime
   t = t + dt;
   step = step + 1;
   fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
      step, convertTo(t - dt, day), convertTo(t, day));

   P0  = value(P_AD); % Previous step pressure
   x0  = value(x_AD);

   % Newton loop
   resNorm = 1e+99;
   iterCounter = 0;
   while (resNorm > tol) && (iterCounter <= maxIter)
      
      %% Define pressure equation
      [rho_AD, rho0] = deal(rhoFunc(P_AD), rhoFunc(P0));
      [PV_AD , PV0]  = deal(PVfunc(P_AD), PVfunc(P0));
      gradz  = grad(G.cells.centroids(:,3));
      v = -(T/mu).*(grad(P_AD) - avg(rho_AD)*g.*gradz);
      pressEq = (1/dt)*(PV_AD.*rho_AD - PV0.*rho0) + div(avg(rho_AD).*v);

      % Adding boundary conditions
      qProd  = -PIprod.*(1/mu).*(P_AD(prodCell) - Pperf);
      % The BC of the production well can be set in 2 different ways, either using the BHP directly or using the equivalent production flowrate.
      pressEq(prodCell) = pressEq(prodCell) - qProd*rho_AD(prodCell); % Method 1
%       pressEq(prodCell) = P_AD(prodCell) - BHPprod;      % Method 2

      pressEq(injCell)  = pressEq(injCell) - qInj*rhoRef;

      %% Define the mass transport equation of the desired species
      % water upstream-index
      upcw = v > 0;
      X_AD = x_AD*MW_Ca/(1*MW_H2O + 0*MW_Ca);
      X0   = x0*MW_Ca/(1*MW_H2O + 0*MW_Ca);
      accumTerm  = (1/dt).*(PV_AD.*rho_AD.*X_AD - PV0.*rho0.*X0);
      % Flux term without dispersion
%       fluxTerm   = div(upwind(X_AD.*rho_AD, upcw).*v);
      % Flux term with dispersion
      D = 10^(-9)*meter^2/second*ones;
      fluxTerm = div(upwind(X_AD.*rho_AD, upcw).*v - avg(rock.poro.*rho_AD).*D.*grad(X_AD));

      transEq = accumTerm + fluxTerm;
      
      %% Add boundary conditions
%       transEq(injCell)  = C_AD(injCell) - influxConc;   % Method 1
      transEq(injCell)  = transEq(injCell) - qInj*rhoRef*Xinflux;   % Method 2
      transEq(prodCell) = transEq(prodCell) - qProd*rho_AD(prodCell)*X_AD(prodCell);

      %% Collect and concatenate all equations (i.e., assemble and linearize system)
      eqs = {pressEq, transEq};
      eq  = cat(eqs{:});
      J   = eq.jac{1};  % Jacobian
      res = eq.val;     % residual
      increment = -(J \ res); % Newton update

      % Update variables
      P_AD.val   = P_AD.val + increment(pIx);
      x_AD.val   = x_AD.val + increment(xIx);

      resNorm = norm(res);
      iterCounter   = iterCounter + 1;
      fprintf('  Iteration %3d:  Res = %.4e\n', iterCounter, resNorm);
   end

   if iterCounter > maxIter
      error('Newton solver did not converge')
   else % Store solution
      sol(step+1)  = struct('time', t, 'pressure', value(P_AD)/psia, 'x', value(x_AD));
      waitbar(t/totTime,hwb);
   end
end
close(hwb);

%% For a 1-dimensional domain
clf
for i = 1:numel(sol)
    plot(1:nc, sol(i).x, 'LineWidth', 1)
    title(['Tracer Concentration Profile (', 'Time Step: ', num2str(i),')'])
    xlim([1,nc])
    ylim([xInit(1), xInflux])
    xlabel('Gridblock Index')
    ylabel('Tracer Concentration [mole fraction]')
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
