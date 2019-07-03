%{
Function to return the forces in a cell (and sum them). Calls
'calcAllForces' while only needing 'cellInfo' as a parameter
('calcAllForces' requires the parameters to be split because it is called
by 'fsolve').
%}
function force = forces(cellInfo)
  force = calcAllForces([cellInfo.xPosition, cellInfo.yPosition],cellInfo);
  sum(abs(force))
end