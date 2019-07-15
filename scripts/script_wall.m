%{
Cell wall scenario script

%}
%% Initialization
dt = 0.0001
dtPlot = dt * 10;
totalTime = 1;

%%
c = initializeNetwork(10, 'wall');
c_def = deformCellDisplacement(c, [], [], [11,12], [11,12]);
plotCell(c_def)

%%
c_new = stepEuler(c_def, dt, totalTime, dtPlot);

%%
plotCell(c_new);