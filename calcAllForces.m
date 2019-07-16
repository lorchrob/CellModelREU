%{
Function to calculate all forces in the network, to be used by 'fsolve'
function. Puts forces into column vector.
%}
function allForces = calcAllForces(positions, cellInfo)
  allForces = zeros(cellInfo.totalNodeCount, 2);
  
  % set positions of fixed nodes
  for i = 1 : numel(positions(:,1))
    if cellInfo.isFixed(i, 1)
      positions(i, 1) = cellInfo.xPosition(i);
    end
    
     if cellInfo.isFixed(i, 2)
      positions(i, 2) = cellInfo.yPosition(i);
    end
  end
  
  % calculate forces of unfixed nodes
  for i = 1 : numel(positions(:,1))
    force = calculateForce(positions, i, cellInfo); 
    if ~cellInfo.isFixed(i, 1)
      allForces(i, 1) = force(1);
    end
   
    if ~cellInfo.isFixed(i, 2)
      allForces(i, 2) = force(2);
    end
  end
end