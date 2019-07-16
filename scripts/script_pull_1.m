%{
Cell is pulled out on all outer nodes with force 'f'
%}
dt = 0.0001;
timeEnd = 2.5;
outerNodeCount = 10;
f = 1000;

dtPlot = dt * 50;

theta = linspace(0,2*pi, outerNodeCount + 1);
theta = theta(1:end-1);
forcesX = cos(theta) * f;
forcesY = sin(theta) * f;

%%
addpath('C:\Users\rlorch\Desktop\CellModelREU')
c = initializeNetwork(outerNodeCount);
c_def = deformCellForce(c, 1:outerNodeCount, [forcesX;forcesY]', [11,12], [11,12]);
plotCell(c_def);

%%
c_new = stepEuler(c_def, dt, timeEnd, dtPlot);

%%
plotCell(c_new)
