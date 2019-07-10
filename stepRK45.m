% Script to initialize a network, deform it, and call ode15s to solve the
% velocity system. Also plots.
%% Init
c = initializeNetwork(10);
c_def = deformCellDisplacement(c, [4], [], [0, 1]);


%% Stepping
pos0 = zeros(c_def.totalNodeCount*2,1);
pos0(1:2:end) = c_def.xPosition;
pos0(2:2:end) = c_def.yPosition;

[t,pos] = ode45(@(t,x) getVelocities(x, c_def), [0 10], pos0);

%% Plotting
 c_new = c;
 c_new.xPosition = pos(end,1:2:end)';
 c_new.yPosition = pos(end,2:2:end)';
 c_new = calculateNodeInfo(c_new);
 plotCell(c_new);