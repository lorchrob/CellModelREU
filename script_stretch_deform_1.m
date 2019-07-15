%{
Script to show cell dyanmics of basic stretched cell.
Cell is displaced by stretching horizontally
%}
dt = 0.001;
timeEnd = 2.5;

dtPlot = dt * 10;

%%
c = initializeNetwork(10);
c_def = deformCellDisplacement(c, [1,6, 2,5, 7,10 , 4,8, 12,14], [ 10,0 ; -10,0 ; 8,-0.5 ; -8,-0.5 ; -8,0.5 ; 8,0.5 ; -6,0; -6,0; -5,0.25 ; -5,-0.25], [3,9], [3,9]);
plotCell(c_def);

%%
c_new = stepEuler(c_def, dt, timeEnd, dtPlot);

%%
plotCell(c_new)

