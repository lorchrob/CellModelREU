% Script to initialize a network, deform it, and call ode15s to solve the
% velocity system. Also plots.
%% Init
c2 = initializeNetwork(10);
c_def2 = deformCellDisplacement(c2, [4], [0,1]);


%% Stepping
pos0 = zeros(c_def2.totalNodeCount*2,1);
pos0(1:2:end) = c_def2.xPosition;
pos0(2:2:end) = c_def2.yPosition;

yp0 = zeros(c_def2.totalNodeCount*2, 1);

%[y0_new, yp0_new] = decic(@)

[t2,pos2] = ode15s(@(t,x) getVelocities(x, c_def2), [0, 7], pos0);

%% Plotting
 c_new2 = c2;
 c_new2.xPosition = pos2(end,1:2:end)';
 c_new2.yPosition = pos2(end,2:2:end)';
 c_new2 = calculateNodeInfo(c_new2);
 plotCell(c_new2);