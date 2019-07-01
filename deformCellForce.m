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
%}
function cellInfoNew = deformCellForce(cellInfo, nodeNums, externalForces)
  cellInfoNew = cellInfo;
  
  for i = 1 : numel(nodeNums)
    % apply external forces to the nodes
    cellInfoNew.externalForces(nodeNums(i),:) = cellInfoNew.externalForces(nodeNums(i),:) + externalForces(i,:);
  end
end