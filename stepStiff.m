% Script to initialize a network, deform it, and call ode15s to solve the
% velocity system. Also plots. 
%% Init
r = initializeNetwork(14);
r_def = deformCellForce(c, 1:25, [2,14], [repmat(-1,25,1), repmat(0,25,1)]);
plotCell(c_def);


%% Stepping
pos0 = zeros(r_def.totalNodeCount*2,1);
pos0(1:2:end) = r_def.xPosition;
pos0(2:2:end) = r_def.yPosition;

yp0 = zeros(r_def.totalNodeCount*2, 1);

%[y0_new, yp0_new] = decic(@)

[t2,pos2] = ode15s(@(t,x) getVelocities(x, r_def), [0, 50], pos0);

%% Plotting
 r_new = r;
 r_new.xPosition = pos2(end,1:2:end)';
 r_new.yPosition = pos2(end,2:2:end)';
 r_new = calculateNodeInfo(r_new);
 plotCell(r_new);