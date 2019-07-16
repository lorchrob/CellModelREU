%{

Cell is trapped in a box

%}
%% Initialization
dt = 0.0001
dtPlot = dt * 10;
totalTime = 0.5;

boxLen = 11;

l = boxLen / 2;

xw = [-l,  l, l, -l, -l];
yw = [-l, -l, l,  l, -l];

%%
addpath('C:\Users\rlorch\Desktop\CellModelREU')
c = initializeNetwork(10, 'wall', xw, yw);
c_def = deformCellDisplacement(c, [], [], [11,12], [11,12]); % 
c_new = stepEuler(c_def, dt, totalTime, dtPlot);

%%
plotCell(c_new);