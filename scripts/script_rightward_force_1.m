%{
Node 1,2,10 have rightward force
%}
dt = 0.001;
timeEnd = 2.5;

dtPlot = dt * 20;

f = 1000; % force

%%
addpath('C:\Users\rlorch\Desktop\CellModelREU')
c = initializeNetwork(10);
c_def = deformCellForce(c, [1,2,10], [f,0 ; f,0 ; f,0], [5,6,7], [5,6,7]);
c.noMeanYChange = true;
plotCell(c_def);

%%
c_new = stepEuler(c_def, dt, timeEnd, dtPlot);

%%
plotCell(c_new)
