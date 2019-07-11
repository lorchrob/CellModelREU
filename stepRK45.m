% Script to initialize a network, deform it, and call ode15s to solve the
% velocity system. Also plots. 
%% Init
externalNodeCount = 15;
r = initializeNetwork(externalNodeCount);
r_def = deformCellForce(r, 1:25, [2,14], [repmat(-100,25,1), repmat(0,25,1)]);


%% Stepping
pos0 = zeros(r_def.totalNodeCount*2,1);
pos0(1:2:end) = r_def.xPosition;
pos0(2:2:end) = r_def.yPosition;
options = odeset('OutputFcn', @(t, y, flag) plotCellWrapper(t, y, flag, externalNodeCount));

[t2,pos2] = ode45(@(t,x) getVelocities(x, r_def), [0, 50], pos0, options);