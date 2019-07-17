%{

%}
%%
dt = 0.00005;
timeEnd = 20;

dtPlot = dt * 20

extNodes = 15;

%%
c = initializeNetwork(extNodes, 'wall'); % 
c_def = deformCellDisplacement(c, [1:27], [repmat(-160, c.totalNodeCount, 1), zeros(c.totalNodeCount,1)], [16, 17], [18], [200, 200], [0] )
c_def.noMeanXChange = false;
c_def.noMeanYChange = false;
plotCell(c_def);

%%
c_new = stepEuler(c_def, dt, timeEnd, dtPlot);

%%
tic;
pos0 = zeros(c_def.totalNodeCount*2,1);
pos0(1:2:end) = c_def.xPosition;
pos0(2:2:end) = c_def.yPosition;
options = odeset('OutputFcn', @(t, y, flag) plotCellWrapper(t, y, flag, extNodes, 'wall', 'timeStepper'));

[t,pos] = ode15s(@(t,x) getVelocities(x, c_def), [0, 1.62], pos0, options);
toc;
%%




