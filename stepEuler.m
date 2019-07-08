%{
Simple function to time step, using Euler's method. t_0 (initial time) is
assumed to be 0).
%}
function cellInfoNew = stepEuler(cellInfo, dt, totalTime) 
  cellInfo.modelType = "timeStepper";
  cellInfoNew = cellInfo;
  tPlot = 0;
  dtPlot = 0.1;
  
  for i = dt:dt:totalTime
%     [A, b] = velSystem(cellInfoNew);
%     vels = A\b;
    pos = zeros(cellInfo.totalNodeCount * 2,1);
    pos(1:2:end) = cellInfoNew.xPosition;
    pos(2:2:end) = cellInfoNew.yPosition;
    vels = getVelocities(pos, cellInfoNew);

    % Euler equations
    cellInfoNew.xPosition = cellInfoNew.xPosition + vels(1:2:end)*dt;
    cellInfoNew.yPosition = cellInfoNew.yPosition + vels(2:2:end)*dt;
    cellInfoNew.xVelocity = vels(1:2:end);
    cellInfoNew.yVelocity = vels(2:2:end);
    cellInfoNew = calculateNodeInfo(cellInfoNew);
    
    % attempt to plot and print at specified time steps
    if i >= tPlot
      plotCell(cellInfoNew);
      drawnow;
      fprintf("Time: %.1f\n", i);
      tPlot = tPlot + dtPlot;
    end
  end
end