%{
Cell and Wall scenario with 20 external nodes (44 total nodes)

node 36 is center node. Makes sense to fix its position
node 35 is top left of center node
node 23 is bottom right of center node

%}
%% Initialization
addpath('C:\Users\rlorch\Desktop\CellModelREU')
dt = 0.0001;
dtPlot = dt * 10;
totalTime = 0.1;

%%
c = initializeNetwork(20, 'wall');
%c_def = deformCellDisplacement(c, [], [], [23, 36], [36]);
c_def = deformCellDisplacement(c, [], [], [36] , [36,1,11]);
c_def.noMeanXChange = true;
%c_def.noMeanYChange = true;
plotCell(c_def)

%%
c_new = stepEuler(c_def, dt, totalTime, dtPlot);


%%
% pos0 = zeros(c_def.totalNodeCount*2,1);
% pos0(1:2:end) = c_def.xPosition;
% pos0(2:2:end) = c_def.yPosition;
% options = odeset('OutputFcn', @(t, y, flag) plotCellWrapper(t, y, flag, 20, 'wall', 'timeStepper'));
% 
% [t,pos] = ode15s(@(t,x) getVelocities(x, c_def), [0, 100], pos0, options);

%%
plotCell(c_new);

%%
% neighbors = [];
% for i = 1:c.totalNodeCount
%   neighbors(i) = numel(c.nodesAdjacent{i})
% end

