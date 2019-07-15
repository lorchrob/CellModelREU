%{
Function to deform a cell by moving around some of the nodes. The nodes to
be moved are specified by 'nodeNums' (an array of integers), and they are each 
moved by the amount specified in 'positionChanges' (an nx2 array of pairs of
floats). 

Update 7/12/2019 : 'fixedNodes' are prescribed a velocity as specified in
'prescVel' argument. If you want those nodes to be fixed, prescribe
velocity of [0,0]

As an example, if the second integer in 'nodeNums' is 5 and the second pair
of floats in 'positionChanges' is [1 -4], then node 5 will be moved 1 to
the right (x direction) and 4 down (y direction).

NOTE: There must be exactly one pair of position changes (x, y) for each node num
(i.e., numel(positionChanges) = 2 * numel(nodeNums)).
%}
function cellInfoNew = deformCellDisplacement(cellInfo, nodeNums, fixedX, fixedY, positionChanges, prescXVel, prescYVel, noMeanXChange, noMeanYChange)
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

  if exist('noMeanXChange', 'var') 
    cellInfo.noMeanXChange = noMeanXChange;
  end
  
  if exist('noMeanYChange', 'var') 
    cellInfo.noMeanYChange = noMeanYChange;
  end
  
%   cellInfoNew.isFixed = false(cellInfo.totalNodeCount, 1);
  cellInfoNew.isFixed(fixedX, 1) = true; 
  cellInfoNew.isFixed(fixedY, 2) = true;
  cellInfoNew.xFix = fixedX;
  cellInfoNew.yFix = fixedY;
end



