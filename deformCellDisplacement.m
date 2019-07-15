%{
Function to deform a cell by moving around some of the nodes. The nodes to
be moved are specified by 'nodeNums' (an array of integers), and they are each 
moved by the amount specified in 'positionChanges' (an nx2 array of pairs of
floats). 

As an example, if the second integer in 'nodeNums' is 5 and the second pair
of floats in 'positionChanges' is [1 -4], then node 5 will be moved 1 to
the right (x direction) and 4 down (y direction).

NOTE: There must be exactly one pair of position changes (x, y) for each node num
(i.e., numel(positionChanges) = 2 * numel(nodeNums)).


Arguments:
  - cellInfo : from running initializeNetwork
  - nodeNums : vector. the nodes that you want to displace
  - positionChanges : the position changes for 'nodeNums'. Matrix (nx2) in the
  form [x1, y1 ; x2, y2 ; x3, y3; ...] where x1 is the change in x position
  of first index of 'nodeNum'.
  - fixedX : the nodes that you want to fix the 'x' position of
  - fixedY : the nodes that you want to fix the 'y' position of
  - prescXVel : the prescribed x velocity for each node specified in
  'fixedX'. 
  - prescYVel : the prescribed y velocity for each node specified in
  'fixedY'.

NOTE:
if you set cellInfo.noMeanXChange to true, then the forces equation for the
FIRST index of fixedX and/or fixedY are replaced by the equation the
specifies that the sum of x(or y) velocities is 0.

%}
function cellInfoNew = deformCellDisplacement(cellInfo, nodeNums, positionChanges, fixedX, fixedY, prescXVel, prescYVel)
  cellInfoNew = cellInfo;

  for i = 1 : numel(nodeNums)
    % displace the nodes
    cellInfoNew.xPosition(nodeNums(i)) = cellInfoNew.xPosition(nodeNums(i)) + positionChanges(i, 1);
    cellInfoNew.yPosition(nodeNums(i)) = cellInfoNew.yPosition(nodeNums(i)) + positionChanges(i, 2);
    
    cellInfoNew = calculateNodeInfo(cellInfoNew);
    % fix (hold constant, not repair) the nodes
    %cellInfoNew.isFixed(nodeNums(i)) = true; 
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



