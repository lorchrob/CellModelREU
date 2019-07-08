%{
Main function for finding the steady state for the system of viscoelastic
elements. Still unsure of parameters/return variables.
%}
function cellInfoNew = findSteadyState(cellInfo) 
  cellInfo.modelType = "steadyStateSolver";
  % set initial guess
  x_0 = [cellInfo.xPosition, cellInfo.yPosition];
  
  options = optimoptions(@fsolve,'MaxFunctionEvaluations', 100000, 'MaxIterations', 10000)
  [newPositions, forceValues] = fsolve(@(x) calcAllForces(x, cellInfo), x_0, options);
  
  % use positions to create cellInfoNew, which is the cell at steady state
  cellInfoNew = cellInfo;
  cellInfoNew.xPosition = newPositions(:,1);
  cellInfoNew.yPosition = newPositions(:,2);
  
  cellInfoNew = calculateNodeInfo(cellInfoNew);
  
  % display force <x, y> at each node
  forceValues
end 



function cross = crossProd(v1, v2)
  cross = v1(1) * v2(2) - v1(2) * v2(1);
end
