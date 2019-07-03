%{
Function to calculate all forces in the network, to be used by 'fsolve'
function. Puts forces into column vector.
%}
function allForces = calcAllForces(positions, cellInfo)
  allForces = zeros(cellInfo.totalNodeCount, 2);
  
  % set positions of fixed nodes
  for i = 1 : numel(positions(:,1))
    if cellInfo.isFixed(i)
      positions(i,:) = [cellInfo.xPosition(i), cellInfo.yPosition(i)];
    end
  end
  
  for i = 1 : numel(positions(:,1))
    % set forces of fixed nodes to 0
    if cellInfo.isFixed(i)
      allForces(i,:) = [0 0];
    % otherwise, calculate force using our equation
    else
      allForces(i,:) = calculateForce(positions, i, cellInfo);
    end
  end
end