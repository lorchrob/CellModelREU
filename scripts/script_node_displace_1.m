%{
1 node is displaced.
%}
tic;
dt = 0.0025;
timeEnd = 2.5;

dtPlot = dt * 5;

%%
c = initializeNetwork(10);
c_def = deformCellDisplacement(c, [4, 9], [0,5 ; 0,-5], [1,6], [1,6]);
plotCell(c_def);

%%
c_new = stepEuler(c_def, dt, timeEnd, dtPlot);

%%
plotCell(c_new)
toc

