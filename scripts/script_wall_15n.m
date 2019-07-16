%{
Cell and Wall scenario with 15 external nodes (27 total nodes)

node 18 is center (top right) node.
node 17 is center (bottom right)
node 16 is center (top left)
node 19 is center (bottom left)

n17: fix x pos
n18: no mean x change
n19: no mean y change

%}
%% Initialization
addpath('C:\Users\rlorch\Desktop\CellModelREU')
dt = 0.00001
dtPlot = dt * 10;
totalTime = 0.5;

%%
c = initializeNetwork(15, 'wall');
c_def = deformCellDisplacement(c, [], [], [18, 17], [19]);
c_def.noMeanXChange = true;
c_def.noMeanYChange = true;
plotCell(c_def)

%% Explicit Solver
c_new = stepEuler(c_def, dt, totalTime, dtPlot);

%% Implicit Solver
pos0 = zeros(c_def.totalNodeCount*2,1);
pos0(1:2:end) = c_def.xPosition;
pos0(2:2:end) = c_def.yPosition;
options = odeset('OutputFcn', @(t, y, flag) plotCellWrapper(t, y, flag, 15, 'wall', 'timeStepper'));

[t,pos] = ode15s(@(t,x) getVelocities(x, c_def), [0, 100], pos0, options);

%%
plotCell(c_new);

%%
% neighbors = [];
% for i = 1:c.totalNodeCount
%   neighbors(i) = numel(c.nodesAdjacent{i})
% end

