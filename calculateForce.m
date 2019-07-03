%{
Function to calculate the total force acting on a node. At steady
state, this function should equal 0.
%}
function force = calculateForce(positions, nodeNum, cellInfoRef)
% first, determine if node is internal or external, and then solve the
% appropriate force equation
  nodePos = [positions(nodeNum, 1), positions(nodeNum, 2)];
  force = [0, 0];
  % force for internal nodes
  if nodeNum > cellInfoRef.externalNodeCount % check if broken
    for i = 1 : numel(cellInfoRef.nodesAdjacent{nodeNum})
      node2 = cellInfoRef.nodesAdjacent{nodeNum}(i);
      %node2Pos = [cellInfoRef.xPosition(node2), cellInfoRef.yPosition(node2)];
      node2Pos = [positions(node2, 1), positions(node2, 2)];
      currentNorm = norm(node2Pos - nodePos);
      force = force + cellInfoRef.k_ti * (currentNorm / cellInfoRef.refLengths{nodeNum}(i) - 1) * (node2Pos - nodePos) / currentNorm;
    end
    
  % force for external nodes 
  else
    for i = 1 : numel(cellInfoRef.nodesAdjacent{nodeNum})
      node2 = cellInfoRef.nodesAdjacent{nodeNum}(i);      
      %node2Pos = [cellInfoRef.xPosition(node2), cellInfoRef.yPosition(node2)];
      node2Pos = [positions(node2, 1), positions(node2, 2)];
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
    nodeIp2 = mod(nodeIp1 + 1, cellInfoRef.externalNodeCount);
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
    %Ip1Pos = [cellInfoRef.xPosition(nodeIp1), cellInfoRef.yPosition(nodeIp1)]; 
    Ip1Pos = [positions(nodeIp1, 1), positions(nodeIp1, 2)];
    %Ip2Pos = [cellInfoRef.xPosition(nodeIp2), cellInfoRef.yPosition(nodeIp2)];
    Ip2Pos = [positions(nodeIp2, 1), positions(nodeIp2, 2)];
    %Im1Pos = [cellInfoRef.xPosition(nodeIm1), cellInfoRef.yPosition(nodeIm1)];
    Im1Pos = [positions(nodeIm1, 1), positions(nodeIm1, 2)];
    %Im2Pos = [cellInfoRef.xPosition(nodeIm2), cellInfoRef.yPosition(nodeIm2)];
    Im2Pos = [positions(nodeIm2, 1), positions(nodeIm2, 2)];
    
    % Getting angles between external elements
    %angleIp1 = acos(dot((nodePos - Ip1Pos),(Ip2Pos - Ip1Pos)) / ( norm(nodePos - Ip1Pos) * norm(Ip2Pos - Ip1Pos) ));
    angleIp1 = atan2( crossProd(nodePos - Ip1Pos, Ip2Pos - Ip1Pos), dot(nodePos - Ip1Pos, Ip2Pos - Ip1Pos));
    %angleI = acos(dot((Im1Pos - nodePos),(Ip1Pos - nodePos)) / ( norm(Im1Pos - nodePos) * norm(Ip1Pos - nodePos) )); % angle at node i
    angleI = atan2( crossProd(Im1Pos - nodePos, Ip1Pos - nodePos), dot(Im1Pos - nodePos, Ip1Pos - nodePos));
    %angleIm1 = acos(dot((Im2Pos - Im1Pos),(nodePos - Im1Pos)) / ( norm(Im2Pos - Im1Pos) * norm(nodePos - Im1Pos) ));
    angleIm1 = atan2( crossProd(Im2Pos - Im1Pos, nodePos - Im1Pos), dot(Im2Pos - Im1Pos, nodePos - Im1Pos));
    
    % normal vector of element I (pointing outwards from cell)
    normI = (Ip1Pos - nodePos)/norm(Ip1Pos - nodePos); 
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
%     force = force + normIm1 * (-2 * cellInfoRef.k_be * (tan(angleI/2) - tan(angleIm1/2))) / (cellInfoRef.refLengths{nodeNum}(end) * norm(nodePos - Im1Pos));
%     force = force + normI * (-2 * cellInfoRef.k_be * (tan(angleIp1/2) - tan(angleI/2))) / (cellInfoRef.refLengths{nodeNum}(1) * norm(Ip1Pos - nodePos));
%     
  end
  
  % add external force, zero vector unless otherwise changed by
  % 'deformCellForce'
  force = force + cellInfoRef.externalForces(nodeNum,:);
end

function cross = crossProd(v1, v2)
  cross = v1(1) * v2(2) - v1(2) * v2(1);
end