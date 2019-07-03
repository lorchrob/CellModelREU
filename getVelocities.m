%{
Using velocity system, returns velocities when given positions
%}
function velocities = getVelocities(positions, cellInfo)
  cellInfoNew = cellInfo;
  cellInfoNew.xPosition = positions(1:2:end);
  cellInfoNew.yPosition = positions(2:2:end);
  cellInfoNew = calculateNodeInfo(cellInfoNew);
  
  [A,b] = velSystem(cellInfoNew);
  velocities = A\b;
end