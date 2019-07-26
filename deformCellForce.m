%{
Function to deform a cell by apply external forces to some of the nodes. The nodes to
apply forces to are specified by 'nodeNums' (an array of integers), and they are each 
pushed/pulled by the amount specified in 'externalForces' (an nx2 array of pairs of
floats). 

Arguments:
  * cellInfo : from running initializeNetwork
  * nodeNums : vector. the nodes that you want to displace
  * externalForces : the forces applied to 'nodeNums', in the form of an
                     nx2 vector (the number of vectors should be the same 
                     as the number of nodeNums)
  * fixedX : the nodes that you want to fix the 'x' position of
  * fixedY : the nodes that you want to fix the 'y' position of
  * prescXVel : the prescribed x velocity for each node specified in
                'fixedX'. 
  * prescYVel : the prescribed y velocity for each node specified in
                'fixedY'.
  * last two arguments: see NOTE below
  * first five arguments REQUIRED, remaining are OPTIONAL
  * must have three pieces of information to solve (essentially you must fix
    at least three x/y positions, ask Barber if confused about this)

NOTE:
if you set cellInfo.noMeanXChange to true, then the forces equation for the
FIRST index of fixedX and/or fixedY are replaced by the equation the
specifies that the sum of x(or y) velocities is 0.
%}
function cellInfoNew = deformCellForce(cellInfo, nodeNums, externalForces, fixedX, fixedY, prescXVel, prescYVel, noMeanXChange, noMeanYChange)

  cellInfoNew = cellInfo;
  
  for i = 1 : numel(nodeNums)
    % apply external forces to the nodes
    cellInfoNew.externalForces(nodeNums(i),:) = cellInfoNew.externalForces(nodeNums(i),:) + externalForces(i,:);
  end
  
  if exist('prescXVel', 'var') && exist('prescYVel', 'var')
    cellInfoNew.xv(fixedX) = prescXVel(:);
    cellInfoNew.yv(fixedY) = prescYVel(:);
  end
  
  if exist('noMeanXChange', 'var') 
    cellInfoNew.noMeanXChange = logical(noMeanXChange);
  end
  
  if exist('noMeanYChange', 'var') 
    cellInfoNew.noMeanYChange = logical(noMeanYChange);
  end
  
%   cellInfoNew.isFixed = false(cellInfo.totalNodeCount, 1);
  cellInfoNew.isFixed(fixedX, 1) = true; 
  cellInfoNew.isFixed(fixedY, 2) = true;
  cellInfoNew.xFix = fixedX;
  cellInfoNew.yFix = fixedY;
  
end
