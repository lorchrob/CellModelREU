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
  if nodeNum > cellInfoRef.externalNodeCount 
    for i = 1 : numel(cellInfoRef.nodesAdjacent{nodeNum})
      node2 = cellInfoRef.nodesAdjacent{nodeNum}(i);
      node2Pos = [positions(node2, 1), positions(node2, 2)];
      currentNorm = norm(node2Pos - nodePos);
      force = force + cellInfoRef.k_ti * (currentNorm / cellInfoRef.refLengths{nodeNum}(i) - 1) * (node2Pos - nodePos) / currentNorm;
    end
    
  % force for external nodes 
  else
    for i = 1 : numel(cellInfoRef.nodesAdjacent{nodeNum})
      node2 = cellInfoRef.nodesAdjacent{nodeNum}(i);      
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
    
    % Find node nums of nearby external nodes (used for angle calculations)
    nodeIp1 = circshift(1:cellInfoRef.externalNodeCount, -1);
    nodeIp1 = nodeIp1(nodeNum);
    nodeIp2 = circshift(1:cellInfoRef.externalNodeCount, -2);
    nodeIp2 = nodeIp2(nodeNum);
    nodeIm1 = circshift(1:cellInfoRef.externalNodeCount, 1);
    nodeIm1 = nodeIm1(nodeNum);
    nodeIm2 = circshift(1:cellInfoRef.externalNodeCount, 2);
    nodeIm2 = nodeIm2(nodeNum);
    
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
    angleIp1 = atan2( crossProd(Ip1Pos - nodePos, Ip2Pos - Ip1Pos), dot(Ip1Pos - nodePos, Ip2Pos - Ip1Pos));
    angleI = atan2( crossProd(nodePos - Im1Pos, Ip1Pos - nodePos), dot(nodePos - Im1Pos, Ip1Pos - nodePos));
    angleIm1 = atan2( crossProd(Im1Pos - Im2Pos, nodePos - Im1Pos), dot(Im1Pos - Im2Pos, nodePos - Im1Pos));
    
    % normal vector of element I (pointing outwards from cell)
    normI = (Ip1Pos - nodePos)/norm(Ip1Pos - nodePos); 
    temp = normI(1);
    normI(1) = normI(2);
    normI(2) = -temp;
    
    normIm1 = (nodePos - Im1Pos) / norm(nodePos - Im1Pos);
    temp = normIm1(1);
    normIm1(1) = normIm1(2);
    normIm1(2) = -temp;
    
    % Adding shearing forces
    % note about l_refs: the external references are always equal (first
    % and last l_ref of any external node)
     shearForce1 = normIm1 * (-2 * cellInfoRef.k_be * (tan(angleI/2) - tan(angleIm1/2))) / (cellInfoRef.refLengths{nodeNum}(end) * norm(nodePos - Im1Pos));
     force = force + shearForce1;
     shearForce2 = -normI * (-2 * cellInfoRef.k_be * (tan(angleIp1/2) - tan(angleI/2))) / (cellInfoRef.refLengths{nodeNum}(1) * norm(Ip1Pos - nodePos));
     force = force + shearForce2;
     
  end
  
  % add external force, zero vector unless otherwise changed by
  % 'deformCellForce' function
  force = force + cellInfoRef.externalForces(nodeNum,:);
  
  % calculate and incorporate the force from the wall
  wallForceMag = 2000; % Magnitude of force from wall
  [dists,dot_dist,overlap,norm_x,norm_y,norm_arc,type] = ...
    dist_from_pt_to_line_segs(nodePos(1), nodePos(2), cellInfoRef.xw, cellInfoRef.yw);
  [mindists,inds] = min(dists,[],1);
  %lininds = sub2ind(size(dists),inds,1:cellInfoRef.totalNodeCount); % fix
  my_eps = 0.1;
  in_or_out = inpolygon(nodePos(1),nodePos(2),cellInfoRef.xw,cellInfoRef.yw);
  ramp_func = @(x,e) (x > -e).*(x < 0).*(x+e).^2.*(e-2*x)./e.^3+(x >= 0);
  fxnwf = (wallForceMag*(norm_x(inds).*...
    ramp_func((1-2*in_or_out').*mindists,my_eps)))';
  fynwf = (wallForceMag*(norm_y(inds).*...
    ramp_func((1-2*in_or_out').*mindists,my_eps)))';
  wallForce = [fxnwf, fynwf];
  force = force + wallForce;
end

function cross = crossProd(v1, v2)
  cross = v1(1) * v2(2) - v1(2) * v2(1);
end