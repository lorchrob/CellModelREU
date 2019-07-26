%{
Using velocity system, returns velocities when given positions. Called by
the ODE solver during time integration.

Takes in:
  - new nodal positions (usually inputed from ODE solver)
  - OLD 'cellInfo' structure (cellInfo from previous timestep)

Outputs:
  - 'cellInfo' structure with updated information
      (i.e. lengths, angles, tangent vectors, normal vectors, external
      forces, etc.)

Note: 
The 'cellInfo' passed in through argument does NOT have updated
positions.
%}
function velocities = getVelocities(positions, cellInfo)
  cellInfoNew = cellInfo;
  cellInfoNew.xPosition = positions(1:2:end);
  cellInfoNew.yPosition = positions(2:2:end);
  cellInfoNew = calculateNodeInfo(cellInfoNew);
  
  cellInfoNew = updateForceCellWall(cellInfoNew, cellInfoNew.xw, cellInfoNew.yw);
  
  [A,b] = velSystem(cellInfoNew);
  velocities = A\b;
end
