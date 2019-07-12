%{
Simple function to time step, using Euler's method. t_0 (initial time) is
assumed to be 0).
%}
function cellInfoNew = stepEuler(cellInfo, dt, totalTime, dtPlot) 
  cellInfo.modelType = "timeStepper";
  cellInfoNew = cellInfo;
  tPlot = 0;
  
  if ~exist('dtPlot', 'var')
    dtPlot = dt * 10;
  end
  
  allPos = zeros(cellInfoNew.totalNodeCount,2);
  step = 1;
  
  for i = dt:dt:totalTime
%     [A, b] = velSystem(cellInfoNew);
%     vels = A\b;
    pos = zeros(cellInfo.totalNodeCount * 2,1);
    pos(1:2:end) = cellInfoNew.xPosition;
    pos(2:2:end) = cellInfoNew.yPosition;
    vels = getVelocities(pos, cellInfoNew);
    
    allPos(:,:,step) = [cellInfoNew.xPosition, cellInfoNew.yPosition];
    step = step + 1;

    % Euler equations
    cellInfoNew.xPosition = cellInfoNew.xPosition + vels(1:2:end)*dt;
    cellInfoNew.yPosition = cellInfoNew.yPosition + vels(2:2:end)*dt;
    cellInfoNew.xVelocity = vels(1:2:end);
    cellInfoNew.yVelocity = vels(2:2:end);
    cellInfoNew = calculateNodeInfo(cellInfoNew);
    
    if norm(vels(1:2:end) * dt) < 0.0001 & norm(vels(2:2:end) * dt) < 0.0001
      break
    end
    
    % attempt to plot and print at specified time steps
    if i >= tPlot
      plotCell(cellInfoNew);
      drawnow;
      fprintf("Time: %.3f\n", i);
      tPlot = tPlot + dtPlot;
    end
  end
end