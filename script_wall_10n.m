%{
Cell and Wall scenario

makes sense to fix x position of internal node 11.
makes (the most) sense to use internal node to specify no mean X and Y
  change

%}
%% Initialization
dt = 0.00005
dtPlot = dt * 10;
totalTime = 0.5;

%%
c = initializeNetwork(10, 'wall');
c_def = deformCellDisplacement(c, [], [], [12,11], [13]);
c_def.noMeanXChange = true;
c_def.noMeanYChange = true;
plotCell(c_def)

%%
c_new = stepEuler(c_def, dt, totalTime, dtPlot);

%%
plotCell(c_new);


