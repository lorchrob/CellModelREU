%{
Main function for finding the steady state for the system of viscoelastic
elements. Still unsure of parameters/return variables.
%}
function cellInfoNew = findSteadyState(cellInfo)  
  % set initial guess
  x_0 = [0 0]
  % set up space for new cell positions
  newPositions = zeros(2, cellInfo.totalNodeCount);
  j = 1;
  
  for i = 1 : cellInfo.totalNodeCount
    newPositions(j:j+1) = fsolve(@(x) calculateForce(x, i, cellInfo), x_0);
    j = j + 2;
  end
  
  % we need to find out how to plot the new state of the cell correctly
  %plot(newPositions(1,:), newPositions(2,:));
  
  % use positions to create cellInfoNew, which is the cell at steady state
  cellInfoNew = cellInfo;
  cellInfoNew.xPosition = transpose(newPositions(1,:));
  cellInfoNew.yPosition = transpose(newPositions(2,:));
  
  % find out how to access 'nodeInfo' function
  cellInfoNew = calculateNodeInfo(cellInfoNew);
end 

%{
Function to calculate the total force acting on a node. At steady
state, this function should equal 0.
%}
function force = calculateForce(nodePos, nodeNum, cellInfoRef)
% first, determine if node is internal or external, and then solve the
% appropriate force equation
  force = [0, 0];
  % force for internal nodes
  if nodeNum > cellInfoRef.externalNodeCount % check if broken
    for i = 1 : numel(cellInfoRef.nodesAdjacent{nodeNum})
      node2 = cellInfoRef.nodesAdjacent{nodeNum}(i);
      node2Pos = [cellInfoRef.xPosition(node2), cellInfoRef.yPosition(node2)];
      currentNorm = norm(node2Pos - nodePos);
      force = force + cellInfoRef.k_ti * (currentNorm / cellInfoRef.refLengths{nodeNum}(i) - 1) * (node2Pos - nodePos) / currentNorm;
    end
  % force for external nodes 
  else
    for i = 1 : numel(cellInfoRef.nodesAdjacent{nodeNum})
      node2 = cellInfoRef.nodesAdjacent{nodeNum}(i);      
      node2Pos = [cellInfoRef.xPosition(node2), cellInfoRef.yPosition(node2)];
      currentNorm = norm(node2Pos - nodePos);
      % first case, edge is external (use k_te)
      if node2 <= cellInfoRef.externalNodeCount
        force = force + cellInfoRef.k_te * (currentNorm / cellInfoRef.refLengths{nodeNum}(i) - 1) * (node2Pos - nodePos) / currentNorm;

      % second case, edge is internal (use k_te)
      else
        force = force + cellInfoRef.k_ti * (currentNorm / cellInfoRef.refLengths{nodeNum}(i) - 1) * (node2Pos - nodePos) / currentNorm;
      end
    end
    
    % Getting node numbers for neighboring external nodes
    % use these nodes in shearing force calculations
    nodeIp1 = mod(nodeNum + 1, cellInfoRef.externalNodeCount);
    if nodeIp1 == 0
      nodeIp1 = cellInfoRef.externalNodeCount;
    end
    nodeIp2 = nodeIp1 + 1;
    if nodeIp2 == 0
      nodeIp2 = cellInfoRef.externalNodeCount;
    end
    nodeIm1 = mod(nodeNum - 1, cellInfoRef.externalNodeCount);
    if nodeIm1 == 0
      nodeIm1 = cellInfoRef.externalNodeCount;
    end
    nodeIm2 = nodeIm1 - 1;
    if nodeIm2 == 0
      nodeIm2 = cellInfoRef.externalNodeCount;
    end
    
    % Getting corresponding node positions for neighboring external nodes
    Ip1Pos = [cellInfoRef.xPosition(nodeIp1), cellInfoRef.yPosition(nodeIp1)]; 
    Ip2Pos = [cellInfoRef.xPosition(nodeIp2), cellInfoRef.yPosition(nodeIp2)];
    Im1Pos = [cellInfoRef.xPosition(nodeIm1), cellInfoRef.yPosition(nodeIm1)];
    Im2Pos = [cellInfoRef.xPosition(nodeIm2), cellInfoRef.yPosition(nodeIm2)];
    
    % Getting angles between external elements
    angleIp1 = acos(dot((nodePos - Ip1Pos),(Ip2Pos - Ip1Pos)) / ( norm(nodePos - Ip1Pos) * norm(Ip2Pos - Ip1Pos) ));
    angleI = acos(dot((Im1Pos - nodePos),(Ip1Pos - nodePos)) / ( norm(Im1Pos - nodePos) * norm(Ip1Pos - nodePos) )); % angle at node i
    angleIm1 = acos(dot((Im2Pos - Im1Pos),(nodePos - Im1Pos)) / ( norm(Im2Pos - Im1Pos) * norm(nodePos - Im1Pos) ));
    
    normI = (Ip1Pos - nodePos)/norm(Ip1Pos - nodePos); % normal vector of element I (pointing outwards from cell)
    temp = normI(1);
    normI(1) = normI(2);
    normI(2) = -temp;
    
    normIm1 = (nodePos - Im1Pos) / norm(nodePos - Im1Pos);
    temp = normIm1(1);
    normIm1(1) = normIm1(2);
    normIm1(2) = -temp;
    
    % Adding shearing stress force
    % note about l_refs: the external references are always equal (first
    % and last l_ref of any external node)
    force = force + normIm1 * (-2 * cellInfoRef.k_be * (tan(angleI/2) - tan(angleIm1/2))) / (cellInfoRef.refLengths{nodeNum}(end) * norm(nodePos - Im1Pos));
    force = force + normI * (-2 * cellInfoRef.k_be * (tan(angleIp1/2) - tan(angleI/2))) / (cellInfoRef.refLengths{nodeNum}(1) * norm(Ip1Pos - nodePos));
    
  end
end
