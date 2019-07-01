%{
Function to calculate more info about the cell (adapted from Barber)
%}
function cellInfo = calculateNodeInfo(cellInfo)
  %  Because this info is stored for each node, some info is redundant. For
  %  instance, each length will get stored twice because each edge has two
  %  nodes.  This is for simplicity and should have negligible impact on
  %  efficiency since the solving process should dominate cpu time.
  cellInfo.lengths = cellInfo.nodesAdjacent;
  cellInfo.dxs = cellInfo.nodesAdjacent;
  cellInfo.dys = cellInfo.nodesAdjacent;
  cellInfo.nxs = cellInfo.nodesAdjacent;
  cellInfo.nys = cellInfo.nodesAdjacent;
 % NOT YET PART OF PROJECT
 % if isfield(cellInfo,'elasticAreas'), firsttime = false; else, firsttime = true; end
  %  These should eventually be estimated in a better way (in particular,
  %  we want a circular shape to correspond to a steady configuration) but
  %  for now we just make them 1s.
  for i = 1:numel(cellInfo.nodesAdjacent)
    cnis = cellInfo.nodesAdjacent{i};
    cellInfo.dxs{i} = cellInfo.xPosition(cnis)-cellInfo.xPosition(i);
    cellInfo.dys{i} = cellInfo.yPosition(cnis)-cellInfo.yPosition(i);
    cellInfo.lengths{i} = sqrt(cellInfo.dxs{i}.^2+cellInfo.dys{i}.^2);
    cellInfo.txs{i} = cellInfo.dxs{i}./cellInfo.lengths{i};
    cellInfo.tys{i} = cellInfo.dys{i}./cellInfo.lengths{i};
    %  These make unit vectors perpendicular to the tangent vectors so that
    %  "normal" crossed with "tangent" yields 1.  The current
    %  counterclockwise orientation of external nodes makes the normal from
    %  external node i to i+1 point outwards (tangent vector points from
    %  node i to i+1)
    cellInfo.nxs{i} = cellInfo.tys{i};
    cellInfo.nys{i} = -cellInfo.txs{i};
    
    % internalRefLength and externalRefLength were used for initial
    % generation, but now for this initialized version, the lengths are the
    % reference lengths (so that there is no tension in the cell)
    
    
%     %  Debugging
%     close(figure(5)); figure(5)
%     quiver(c.xn(i)*ones(size(cnis)),c.yn(i)*ones(size(cnis)),...
%       c.txs{i}'.*c.lengths{i}',c.tys{i}'.*c.lengths{i}',0);
%     hold on
%     midxs = (c.xn(i)*ones(size(cnis))+c.xn(cnis)')./2;
%     midys = (c.yn(i)*ones(size(cnis))+c.yn(cnis)')./2;
%     plot(c.xn(cnis),c.yn(cnis),'x');
%     quiver(midxs,midys,c.nxs{i}',c.nys{i}');
%     axis equal;
      
    % Between each pair of connections, there is an angle, which we
    % calculate:
    cnisp1 = circshift(1:numel(cellInfo.txs{i}),-1); % apply a shift of elements to easily compare this node with the next node
    cellInfo.crosses{i} = cellInfo.txs{i}.*cellInfo.tys{i}(cnisp1)-...
      cellInfo.tys{i}.*cellInfo.txs{i}(cnisp1);
    cellInfo.dots{i} = cellInfo.txs{i}.*cellInfo.txs{i}(cnisp1)+...
      cellInfo.tys{i}.*cellInfo.tys{i}(cnisp1);
    cellInfo.alphs{i} = atan2(cellInfo.crosses{i},cellInfo.dots{i});
    cellInfo.alphs{i} = cellInfo.alphs{i}+2*pi*(cellInfo.alphs{i}<0);
    sum(cellInfo.alphs{i});
    %  Between each pair of connections, there is also an element with
    %  associated info (like area of that element)
    if i <= cellInfo.externalNodeCount
      axs = cellInfo.xPosition(cnis(1:end-1),:); ays = cellInfo.yPosition(cnis(1:end-1),:);
      bxs = cellInfo.xPosition(cnis(2:end),:); bys = cellInfo.yPosition(cnis(2:end),:);
    else
      axs = cellInfo.xPosition(cnis,:); ays = cellInfo.yPosition(cnis,:);
      bxs = cellInfo.xPosition(circshift(cnis,-1),:); bys = cellInfo.yPosition(circshift(cnis,-1),:);
    end
    xxi = cellInfo.xPosition(i)*ones(size(axs)); xyi = cellInfo.yPosition(i)*ones(size(bxs));
    a = [axs,ays]; b = [bxs,bys]; xi = [xxi,xyi];
    
   % NOT YET PART OF PROJECT
   % if firsttime
      % cellInfo.elasticAreas{i} = triangleAreaInfo(xi,a,b);
   % else
      %  c.anfs-"normalized forces".  This is the "geometric portion" of
      %  the area-related forces and multiplying by ka, the area elastic
      %  force modulus, will recover the actual forces
      % [cellInfo.elasticAreas{i},cellInfo.anfs{i}] = ...
      %  triareainfo(xi,a,b,cellInfo.earearefs{i},1);
   % end
  end
  %  Store area and "volume estimate"
  cellInfo.area = polyarea(cellInfo.xPosition(1:cellInfo.externalNodeCount),cellInfo.yPosition(1:cellInfo.externalNodeCount));
  
  if ~isfield(cellInfo, "refLengths")
    cellInfo.refLengths = cellInfo.lengths;
  end
end