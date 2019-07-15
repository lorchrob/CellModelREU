% Script to initialize a network, deform it, and call ode15s to solve the
% velocity system. Also plots. 
%% Init
externalNodeCount = 9;
r = initializeNetwork(externalNodeCount);
r_def = deformCellDisplacement(r, [1 5], [1 5], [1 5], [3 3; -2 -1]);


%% Stepping
pos0 = zeros(r_def.totalNodeCount*2,1);
pos0(1:2:end) = r_def.xPosition;
pos0(2:2:end) = r_def.yPosition;
options = odeset('OutputFcn', @(t, y, flag) plotCellWrapper(t, y, flag, externalNodeCount, 'noWall', 'timeStepper'));

[t,pos] = ode15s(@(t,x) getVelocities(x, r_def), [0, 100], pos0, options);

%% Storage
 r_new = r;
 r_new.xPosition = pos(end,1:2:end)';
 r_new.yPosition = pos(end,2:2:end)';
 r_new = calculateNodeInfo(r_new);