%{
The system of equations for velocity (given position). Returns A, b, where velocity = x and
we have Ax = b. Used by various stepping methods.
%}
function [A,b] = velSystem(cellInfo)
  matSize = cellInfo.totalNodeCount * 2;
  A = zeros(matSize);
  b = zeros(matSize, 1);
  extLineSegs = cellInfo.externalLineSegments; % all the external line segments
  
  
  % Setting the diagonal elements (u1, v1, u2, v2, u3, v3, ...)
  for i = 1:matSize
    node = ceil(i/2);
    adjNodes = cellInfo.nodesAdjacent{node};
    
    % creates a vector of line segments around node 
    lineSegs = [ repmat( node,1,numel(adjNodes) ) ; adjNodes]';
    
    % Vector of booleans saying if each line segment (from above) are
    %   external (1) or internal (0)
    IorE = ismember(lineSegs, extLineSegs, 'rows') + ismember(flip(lineSegs,2), extLineSegs, 'rows');

    XorY = mod(i,2) == 1; % true if 'i' is odd. false if even.
    
    A(i,i) = A(i,i) - sum( ( cellInfo.mu_e * IorE + cellInfo.mu_i * ~IorE ) ./ cellInfo.lengths{node} .* (cellInfo.txs{node}.^2 * XorY + cellInfo.tys{node}.^2 * ~XorY) );
    
    % experimenting with trying to vectorize 'i' loops. Probably impossible
    %   to do since cellInfo is a cell and you can't do vector stuff with
    %   cells.
%   index = 1:2*(matSize+1):matSize^2; % index of odd diagonal elements
%   [I,J] = ind2sub([matSize, matSize], index); % i,j indexes of odd diagonal elements
%    IorE = ismember(cellInfo.externalLineSegments, [
%   A(index) = A(index) - ( ) ./ cellInfo.lengths{ ceil(I/2) } .* cellInfo.txs{ ceil(I/2) } ^ 2 
    
  end
  
  % Setting (some of) the subdiagonal and superdiagonal elements
  for i = 2:2:matSize
    node = ceil(i/2);
    adjNodes = cellInfo.nodesAdjacent{node};
    
    lineSegs = [ repmat(node,1,numel(adjNodes)) ; adjNodes ]';
    IorE = ismember(lineSegs, extLineSegs, 'rows') + ismember(flip(lineSegs,2), extLineSegs, 'rows');
    A(i, i-1) = A(i, i-1) - sum(( cellInfo.mu_e * IorE + cellInfo.mu_i * ~IorE ) ./ cellInfo.lengths{node} .* cellInfo.tys{node} .* cellInfo.txs{node} );
    A(i-1, i) = A(i, i-1);
  end
  
  % Setting elements for (odd rows and odd columns) and (even rows and even columns)
  for i = 1:2:matSize
    node = ceil(i/2);
    adjNodes = cellInfo.nodesAdjacent{node};
    
    lineSegs = [repmat(node,1,numel(adjNodes)) ; adjNodes]';
    IorE = ismember(lineSegs, extLineSegs, 'rows') + ismember(flip(lineSegs,2), extLineSegs, 'rows');

    A(i, adjNodes*2-1) = A(i, adjNodes*2-1) + (( cellInfo.mu_e * IorE + cellInfo.mu_i * ~IorE ) ./ cellInfo.lengths{node} .* cellInfo.txs{node} .^ 2)';
    A(i+1, adjNodes*2) = A(i+1, adjNodes*2) + (( cellInfo.mu_e * IorE + cellInfo.mu_i * ~IorE ) ./ cellInfo.lengths{node} .* cellInfo.tys{node} .^ 2)'; 
  end 
  
  % Setting elements for everything else in A
  for i = 1:2:matSize
    node = ceil(i/2);
    adjNodes = cellInfo.nodesAdjacent{node};
    
    lineSegs = [repmat(node,1,numel(adjNodes)) ; adjNodes]';
    IorE = ismember(lineSegs, extLineSegs, 'rows') + ismember(flip(lineSegs,2), extLineSegs, 'rows');
    A(i+1, adjNodes*2-1) = A(i+1, adjNodes*2-1) + (( cellInfo.mu_e * IorE + cellInfo.mu_i * ~IorE ) ./ cellInfo.lengths{node} .* cellInfo.txs{node} .* cellInfo.tys{node})';
    A(i, adjNodes*2) = A(i+1, adjNodes*2-1);
  end  
    
  % Setting b vector.
  for i = 1:2:matSize
    node = ceil(i/2);
    adjNodes = cellInfo.nodesAdjacent{node};
    
    lineSegs = [repmat(node,1,numel(adjNodes)) ; adjNodes]';
    IorE = ismember(lineSegs, extLineSegs, 'rows') + ismember(flip(lineSegs,2), extLineSegs, 'rows');
    
    tensionX = - sum(( cellInfo.k_te * IorE + cellInfo.k_ti * ~IorE ) .* (cellInfo.lengths{node} ./ cellInfo.refLengths{node} - 1) .* cellInfo.txs{node});
    tensionY = - sum(( cellInfo.k_te * IorE + cellInfo.k_ti * ~IorE ) .* (cellInfo.lengths{node} ./ cellInfo.refLengths{node} - 1) .* cellInfo.tys{node});
    
    b(i) = b(i) + tensionX;
    b(i+1) = b(i+1) + tensionY;
    
    % shearing forces    
    if node <= cellInfo.externalNodeCount
      shiftNodep1 = circshift(1:cellInfo.externalNodeCount, -1);
      shiftNodem1 = circshift(1:cellInfo.externalNodeCount, 1);
      shearForceX = - (-2) * cellInfo.k_be * ( tan((cellInfo.alphs{node}(end)-pi)/2) - tan((cellInfo.alphs{shiftNodem1(node)}(end)-pi)/2) ) / ( cellInfo.refLengths{node}(end) * cellInfo.lengths{node}(end) ) * -cellInfo.nxs{node}(end)...
                    - (-2) * cellInfo.k_be * ( tan((cellInfo.alphs{shiftNodep1(node)}(end)-pi)/2) - tan((cellInfo.alphs{node}(end)-pi)/2) ) / ( cellInfo.refLengths{node}(1) * cellInfo.lengths{node}(1) ) * -cellInfo.nxs{node}(1);
      shearForceY = - (-2) * cellInfo.k_be * ( tan((cellInfo.alphs{node}(end)-pi)/2) - tan((cellInfo.alphs{shiftNodem1(node)}(end)-pi)/2) ) / ( cellInfo.refLengths{node}(end) * cellInfo.lengths{node}(end) ) * -cellInfo.nys{node}(end)...
                    - (-2) * cellInfo.k_be * ( tan((cellInfo.alphs{shiftNodep1(node)}(end)-pi)/2) - tan((cellInfo.alphs{node}(end)-pi)/2) ) / ( cellInfo.refLengths{node}(1) * cellInfo.lengths{node}(1) ) * -cellInfo.nys{node}(1);
      b(i) = b(i) + shearForceX; %+ shearForceX; % also had - for both at one point
      b(i+1) = b(i+1) + shearForceY; %+ shearForceY;
    end
  end
  
  % prescibe external force
  b(1:2:end) = b(1:2:end) - cellInfo.externalForces(:,1) - cellInfo.xwf;
  b(2:2:end) = b(2:2:end) - cellInfo.externalForces(:,2) - cellInfo.ywf;    

  % Fixing nodes that were specified in 'deform...' functions
  A(1:2:end,:) = A(1:2:end,:) .* ~cellInfo.isFixed;
  A(2:2:end,:) = A(2:2:end,:) .* ~cellInfo.isFixed;
  A(1:2:end, 1:2:end) = A(1:2:end, 1:2:end) + eye(cellInfo.totalNodeCount) .* cellInfo.isFixed;
  A(2:2:end, 2:2:end) = A(2:2:end, 2:2:end) + eye(cellInfo.totalNodeCount) .* cellInfo.isFixed;

  b(1:2:end) = b(1:2:end) .* ~cellInfo.isFixed + cellInfo.xv;
  b(2:2:end) = b(2:2:end) .* ~cellInfo.isFixed + cellInfo.yv;

end


