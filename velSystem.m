%{
The system of equations for velocity (given position). Returns A, b, where velocity = x and
we have Ax = b. Used by various stepping methods.
%}
function [A,b] = velSystem(cellInfo)
  matSize = cellInfo.totalNodeCount * 2;
  A = zeros(matSize);
  b = zeros(matSize, 1);
  
  % Setting the odd diagonal elements (u1, u2, u3,...)
  for i = 1:2:matSize
    node = ceil(i/2);
    for j = 1:numel(cellInfo.nodesAdjacent{node})

      % Using mu_external if external element
      if cellInfo.nodesAdjacent{node}(j) <= cellInfo.externalNodeCount & node <= cellInfo.externalNodeCount
        A(i,i) = A(i,i) - cellInfo.mu_e / cellInfo.lengths{node}(j) * cellInfo.txs{node}(j) ^ 2;  
      % Using mu_internal if internal element
      else
        A(i,i) = A(i,i) - cellInfo.mu_i / cellInfo.lengths{node}(j) * cellInfo.txs{node}(j) ^ 2;
      end

    end
  end
  
  % Setting the even diagonal elements (v1, v2, v3,...)
  for i = 2:2:matSize
    node = ceil(i/2);
    for j = 1:numel(cellInfo.nodesAdjacent{node})

      % Using mu_external if external element
      if cellInfo.nodesAdjacent{node}(j) <= cellInfo.externalNodeCount & node <= cellInfo.externalNodeCount
        A(i,i) = A(i,i) - cellInfo.mu_e / cellInfo.lengths{node}(j) * cellInfo.tys{node}(j) ^ 2;  
      % Using mu_internal if internal element
      else
        A(i,i) = A(i,i) - cellInfo.mu_i / cellInfo.lengths{node}(j) * cellInfo.tys{node}(j) ^ 2;
      end

    end
  end
  
  % Setting the subdiagonal and superdiagonal elements
  for i = 2:matSize-1
    node = ceil(i/2);
    for j = 1:numel(cellInfo.nodesAdjacent{node})
      
      % Using mu_external if external element
      if cellInfo.nodesAdjacent{node}(j) <= cellInfo.externalNodeCount & node <= cellInfo.externalNodeCount
        A(i,i-1) = A(i,i-1) - cellInfo.mu_e / cellInfo.lengths{node}(j) * cellInfo.tys{node}(j) * cellInfo.txs{node}(j);
        A(i-1,i) = A(i,i-1); % sub and super diagonals are equal
      % Using mu_internal if internal element
      else
        A(i,i-1) = A(i,i-1) - cellInfo.mu_i / cellInfo.lengths{node}(j) * cellInfo.tys{node}(j) * cellInfo.txs{node}(j);
        A(i-1,i) = A(i,i-1); % sub and super diagonals are equal
      end
     
    end
  end
  
  % Setting elements for odd rows and odd columns
  for i = 1:2:matSize
    node = ceil(i/2);
    for j = 1:numel(cellInfo.nodesAdjacent{node}) % looping through all connected elements
      adjacentNode = cellInfo.nodesAdjacent{node}(j);
      
      % Using mu_external if external element
      if cellInfo.nodesAdjacent{node}(j) <= cellInfo.externalNodeCount & node <= cellInfo.externalNodeCount
        A(i,adjacentNode*2-1) = A(i,adjacentNode*2-1) + cellInfo.mu_e / cellInfo.lengths{node}(j) * cellInfo.txs{node}(j)^2;
      % Using mu_internal if internal element
      else
        A(i,adjacentNode*2-1) = A(i,adjacentNode*2-1) + cellInfo.mu_i / cellInfo.lengths{node}(j) * cellInfo.txs{node}(j)^2;
      end
      
    end
  end
  
  % Setting elements for even rows and even columns
  for i = 2:2:matSize
    node = ceil(i/2);
    for j = 1:numel(cellInfo.nodesAdjacent{node}) % looping through all connected elements
      adjacentNode = cellInfo.nodesAdjacent{node}(j);
      
      % Using mu_external if external element
      if cellInfo.nodesAdjacent{node}(j) <= cellInfo.externalNodeCount & node <= cellInfo.externalNodeCount
        A(i,adjacentNode*2) = A(i,adjacentNode*2) + cellInfo.mu_e / cellInfo.lengths{node}(j) * cellInfo.tys{node}(j)^2;
      % Using mu_internal if internal element
      else
        A(i,adjacentNode*2) = A(i,adjacentNode*2) + cellInfo.mu_i / cellInfo.lengths{node}(j) * cellInfo.tys{node}(j)^2;
      end
      
    end
  end
  
  % Setting elements for
  for i = 1:2:matSize
    node = ceil(i/2);
    for j = 1:numel(cellInfo.nodesAdjacent{node})
      adjacentNode = cellInfo.nodesAdjacent{node}(j);
      
      if cellInfo.nodesAdjacent{node}(j) <= cellInfo.externalNodeCount & node <= cellInfo.externalNodeCount
        A(i+1, adjacentNode*2-1) = A(i+1, adjacentNode*2-1) + cellInfo.mu_e / cellInfo.lengths{node}(j) * cellInfo.txs{node}(j) * cellInfo.tys{node}(j);
        A(i, adjacentNode*2) = A(i+1, adjacentNode*2-1);
      else
        A(i+1, adjacentNode*2-1) = A(i+1, adjacentNode*2-1) + cellInfo.mu_i / cellInfo.lengths{node}(j) * cellInfo.txs{node}(j) * cellInfo.tys{node}(j);
        A(i, adjacentNode*2) = A(i+1, adjacentNode*2-1);
      end
      
    end
  end
  
%   cellInfo.shears = zeros(cellInfo.externalNodeCount, 2);
%   cellInfo.tensions = cellInfo.lengths;
    
  % Setting b vector.
  for i = 1:2:matSize
    node = ceil(i/2);
    for j = 1:numel(cellInfo.nodesAdjacent{node})
      
      if cellInfo.nodesAdjacent{node}(j) <= cellInfo.externalNodeCount & node <= cellInfo.externalNodeCount
        tensionX = - cellInfo.k_te * (cellInfo.lengths{node}(j) / cellInfo.refLengths{node}(j) - 1) * cellInfo.txs{node}(j);
        b(i) = b(i) + tensionX;
        tensionY = - cellInfo.k_te * (cellInfo.lengths{node}(j) / cellInfo.refLengths{node}(j) - 1) * cellInfo.tys{node}(j);
        b(i+1) = b(i+1) + tensionY;
        
%         cellInfo.tensions{node}(j) = [tensionX, tensionY];
        
      else
        tensionX = - cellInfo.k_ti * (cellInfo.lengths{node}(j) / cellInfo.refLengths{node}(j) - 1) * cellInfo.txs{node}(j);
        b(i) = b(i) + tensionX;
        tensionY = - cellInfo.k_ti * (cellInfo.lengths{node}(j) / cellInfo.refLengths{node}(j) - 1) * cellInfo.tys{node}(j);
        b(i+1) = b(i+1) + tensionY;
      end
      
    end
    
    %shearing forces
    
     if node <= cellInfo.externalNodeCount
       shiftNodep1 = circshift(1:cellInfo.externalNodeCount, -1);
       shiftNodem1 = circshift(1:cellInfo.externalNodeCount, 1);
       shearForceX = - (-2) * cellInfo.k_be * ( tan((cellInfo.alphs{node}(end)-pi)/2) - tan((cellInfo.alphs{shiftNodem1(node)}(end)-pi)/2) ) / ( cellInfo.refLengths{node}(1) * cellInfo.lengths{node}(1) ) * -cellInfo.nxs{node}(1)...
                     - (-2) * cellInfo.k_be * ( tan((cellInfo.alphs{shiftNodep1(node)}(end)-pi)/2) - tan((cellInfo.alphs{node}(end)-pi)/2) ) / ( cellInfo.refLengths{node}(end) * cellInfo.lengths{node}(end) ) * -cellInfo.nxs{node}(end);
       shearForceY = - (-2) * cellInfo.k_be * ( tan((cellInfo.alphs{node}(end)-pi)/2) - tan((cellInfo.alphs{shiftNodem1(node)}(end)-pi)/2) ) / ( cellInfo.refLengths{node}(1) * cellInfo.lengths{node}(1) ) * -cellInfo.nys{node}(1)...
                     - (-2) * cellInfo.k_be * ( tan((cellInfo.alphs{shiftNodep1(node)}(end)-pi)/2) - tan((cellInfo.alphs{node}(end)-pi)/2) ) / ( cellInfo.refLengths{node}(end) * cellInfo.lengths{node}(end) ) * -cellInfo.nys{node}(end);
       b(i) = b(i) + shearForceX; %+ shearForceX; % also had - for both at one point
       b(i+1) = b(i+1) - shearForceY; %+ shearForceY;
     end
  end
  
  for i = 1:cellInfo.totalNodeCount
    if cellInfo.isFixed(i)
      A(i*2-1,:) = 0;
      A(i*2-1, i*2-1) = 1;
      
      A(i*2,:) = 0;
      A(i*2, i*2) = 1;
      
      b(i*2-1) = 0;
      b(i*2) = 0;
    
    end
  end
  
end