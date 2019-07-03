%{
Simple function to time step, using Euler's method. t_0 (initial time) is
assumed to be 0).
%}
function cellInfoNew = EulerStep(cellInfo, dt, totalTime) 
  cellInfoNew = cellInfo;
  tPlot = 0;
  dtPlot = 0.5;
  
  for i = dt:dt:totalTime
%     [A, b] = velSystem(cellInfoNew);
%     vels = A\b;
    pos = zeros(cellInfo.totalNodeCount * 2,1);
    pos(1:2:end) = cellInfoNew.xPosition;
    pos(2:2:end) = cellInfoNew.yPosition;
    vels = getVelocities(pos, cellInfoNew);

    cellInfoNew.xPosition = cellInfoNew.xPosition + vels(1:2:end)*dt;
    cellInfoNew.yPosition = cellInfoNew.yPosition + vels(2:2:end)*dt;
%     cellInfoNew = calculateNodeInfo(cellInfoNew);
    
    if i >= tPlot
      plot = cellInfoNew;
      plotCell(plot);
      fprintf("Time: %.1f\n", i);
      tPlot = tPlot + dtPlot;
    end
  end
end