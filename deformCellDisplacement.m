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

NOTE: The displaced nodes will be fixed at their new position.
%}
function cellInfoNew = deformCellDisplacement(cellInfo, nodeNums, positionChanges)
  cellInfoNew = cellInfo;

  for i = 1 : numel(nodeNums)
    % displace the nodes
    cellInfoNew.xPosition(nodeNums(i)) = cellInfoNew.xPosition(nodeNums(i)) + positionChanges(i, 1);
    cellInfoNew.yPosition(nodeNums(i)) = cellInfoNew.yPosition(nodeNums(i)) + positionChanges(i, 2);
    
    cellInfoNew = calculateNodeInfo(cellInfoNew);
    % fix (hold constant, not repair) the nodes
    %cellInfoNew.isFixed(nodeNums(i)) = true; 
  end
  
  % fixing for tests
%   cellInfoNew.isFixed = ones(cellInfo.totalNodeCount, 1);
%   cellInfoNew.isFixed(nodeNums) = false;
  
%   cellInfoNew.isFixed([1,2,3, 9,10]) = true;
  cellInfoNew.isFixed([6]) = true;
end

% c_def = deformCellDisplacement(c, [1,4], [0, 3; 3 -3;]);
% c_def = deformCellDisplacement(c, [1,4, 11], [0, 3; 3 -3; -2 -2]);

