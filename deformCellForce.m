%{
Function to deform a cell by apply external forces to some of the nodes. The nodes to
apply forces to are specified by 'nodeNums' (an array of integers), and they are each 
pushed/pulled by the amount specified in 'externalForces' (an nx2 array of pairs of
floats). 

As an example, if the second integer in 'nodeNums' is 5 and the second pair
of floats in 'externalForces' is [1 -4], then the a force of 1 in the x
direction and -4 in the y direction will be applied. 

NOTE: There must be exactly one pair of forces (x, y) for each node num
(i.e., numel(externalForces) = 2 * numel(nodeNums)).

Arguments:
  - cellInfo : from running initializeNetwork
  - nodeNums : vector. the nodes that you want to displace
  - fixedX : the nodes that you want to fix the 'x' position of
  - fixedY : the nodes that you want to fix the 'y' position of
  - positionChanges : the position changes for 'nodeNums'. Matrix (nx2) in the
  form [x1, y1 ; x2, y2 ; x3, y3; ...] where x1 is the change in x position
  of first index of 'nodeNum'.
  - prescXVel : the prescribed x velocity for each node specified in
  'fixedX'. 
  - prescYVel : the prescribed y velocity for each node specified in
  'fixedY'.

NOTE:
if you set cellInfo.noMeanXChange to true, then the forces equation for the
FIRST index of fixedX and/or fixedY are replaced by the equation the
specifies that the sum of x(or y) velocities is 0.

%}
function cellInfoNew = deformCellForce(cellInfo, nodeNums, externalForces, fixedX, fixedY, prescXVel, prescYVel)
  cellInfoNew = cellInfo;
  
  for i = 1 : numel(nodeNums)
    % apply external forces to the nodes
    cellInfoNew.externalForces(nodeNums(i),:) = cellInfoNew.externalForces(nodeNums(i),:) + externalForces(i,:);
  end
  
  if exist('prescXVel', 'var') && exist('prescYVel', 'var')
    cellInfoNew.xv(fixedX) = prescXVel(:);
    cellInfoNew.yv(fixedY) = prescYVel(:);
  end
  
%   cellInfoNew.isFixed = false(cellInfo.totalNodeCount, 1);
  cellInfoNew.isFixed(fixedX, 1) = true; 
  cellInfoNew.isFixed(fixedY, 2) = true;
  cellInfoNew.xFix = fixedX;
  cellInfoNew.yFix = fixedY;
  
end
