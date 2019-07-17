% Script to initialize a network, deform it, and call ode15s to solve the
% velocity system. Also plots. 
%% Init 1 (squeezed by walls on top and bottom)
externalNodeCount = 20;
r = initializeNetwork(externalNodeCount, 'wall');
r_def = deformCellDisplacement(r, [], [], [21], [21 1 11]);

%% Init 2 (stretched between 2 fixed points)
externalNodeCount = 10;
r = initializeNetwork(externalNodeCount);
r_def = deformCellDisplacement(r, [1 5], [4 4; -3.5 4.5], [1 5], [1 5]);

%% Init 3 (push to one side)
externalNodeCount = 10;
r = initializeNetwork(externalNodeCount);
r_def = deformCellForce(r, 1:14, repmat([-100 0], 14, 1), [1 10], [1 10]);

%% Init 4 
externalNodeCount = 20;
r = initializeNetwork(externalNodeCount, 'wall');
r_def = deformCellDisplacement(r, [], [], [21], [21 1 11]);

%% Init 5 (simple deformation)
externalNodeCount = 10;
r = initializeNetwork(externalNodeCount);
r_def = deformCellDisplacement(r, [2], [3 3], [1 4], [1 4]);

%% Init 6 (deformation with internal node moved)
externalNodeCount = 10;
r = initializeNetwork(externalNodeCount);
r_def = deformCellDisplacement(r, [2], [3 3], [1 4], [1 4]);

%% Stepping
pos0 = zeros(r_def.totalNodeCount*2,1);
pos0(1:2:end) = r_def.xPosition;
pos0(2:2:end) = r_def.yPosition;
options = odeset('OutputFcn', @(t, y, flag) plotCellWrapper(t, y, flag, externalNodeCount, 'wall', 'timeStepper'));

[t,pos] = ode15s(@(t,x) getVelocities(x, r_def), [0, 100], pos0, options);

%% Storage
 r_new = r;
 r_new.xPosition = pos(end,1:2:end)';
 r_new.yPosition = pos(end,2:2:end)';
 r_new = calculateNodeInfo(r_new);