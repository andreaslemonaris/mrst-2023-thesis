%% Define grid and petrophysical properties
G = cartGrid([10, 10, 1], [200, 200, 10]);
G = computeGeometry(G);
rock = struct('perm', ones(G.cells.num, 1) * 100 * milli * darcy, ...
              'poro', ones(G.cells.num, 1) * 0.3);

%% Initialize state
state = initResSol(G, 100 * barsa, [0.2, 0.8]);

% Define fluid properties
fluid = initSimpleADIFluid('phases', 'WOG', 'mu', [1, 10, 0.1] * centi * poise, ...
                           'rho', [1000, 700, 500] * kilogram / meter^3, ...
                           'n', [2, 2, 2]);

% Add polymer component
% fluid = addComponent(fluid, 'polymer', 'concentration', 0);

% Define wells and schedule
W = verticalWell([], G, rock, 5, 5, 1, 'Type', 'rate', 'Val', 1*meter^3/day, ...
                 'Radius', 0.1, 'Name', 'Injector');
schedule = simpleSchedule(repmat(1*day, 10, 1), 'W', W);

%% Run simulation
model = ThreePhaseBlackOilPolymerModel(G, rock, fluid);
[wellSols, states] = simulateScheduleAD(state, model, schedule);

%% Extract and visualize results
polymer_concentration = arrayfun(@(x) mean(x.c), states);
time = cumsum(schedule.step.val);
plot(time, polymer_concentration);
xlabel('Time [days]');
ylabel('Average Polymer Concentration');
