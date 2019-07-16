%{
Prescribed velocity on 3 rightmost nodes
%}
dt = 0.01;
timeEnd = 5;

dtPlot = dt * 10;
v = 5;

%%
addpath('C:\Users\rlorch\Desktop\CellModelREU')
c = initializeNetwork(10);
c_def = deformCellDisplacement(c, [],[], [1,2,10], [1,2,10], [v,v,v],[0,0,0]);
plotCell(c_def);

%%
c_new = stepEuler(c_def, dt, timeEnd, dtPlot);

%%
plotCell(c_new)
