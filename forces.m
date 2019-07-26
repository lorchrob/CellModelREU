%{
Function to return the forces in a cell (and sum them). Calls
'calcAllForces' while only needing 'cellInfo' as a parameter
('calcAllForces' requires the parameters to be split because it is called
by 'fsolve').

Essentially, this function makes it easier to see the forces in a cell if
the user wants to for whatever reason.
%}
function force = forces(cellInfo)
  force = calcAllForces([cellInfo.xPosition, cellInfo.yPosition],cellInfo);
  sum(abs(force))
end