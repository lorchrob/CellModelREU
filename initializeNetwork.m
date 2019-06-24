%{ 
TO DO:
  * make structs to store information about the network
      - number of internal and external nodes
      - position   
      - physical properties such as reference lengths, pressure,
        viscoelasticity, etc.
  * generate mesh
  * find out what is given/passed into the system besides node count
  * possibly use varargin for optional parameters
%}

%{
Main function to initialize the system
NOTE: not sure yet if necessary to use all three returned structs (as in
Barber's program)
NOTE: if I understand correctly, the # of external nodes is specified,
outline is made, and then middle is filled in?
%}
function cellInfo = initializeNetwork(externalNodeCount)  
  % STUB
  % set initial values for all fields of the 'cellInfo' struct
  % NOTE: still unsure which values are given and which are calculated
  s = "not yet implemented"
  cellInfo.externalNodeCount = externalNodeCount;
  cellInfo.radius = 12.4/2; 
  cellInfo.refArea = pi*cellInfo.radius^2;
  
  % distance from center to the external nodes
  cellInfo.internalRefLength = sqrt(2*pi*cellInfo.radius^2/cellInfo.externalNodeCount/sin(2*pi/cellInfo.externalNodeCount)); 
  cellInfo.externalRefLength = 2*sin(2*pi/cellInfo.externalNodeCount/2)*cellInfo.internalRefLength;
  
  % set up outer circle
  nodeAngles = linspace(0, 2*pi, cellInfo.externalNodeCount+1)';
  nodeAngles = nodeAngles(1:end-1);
  cellInfo.xPosition = bsxfun(@plus, 0, cellInfo.internalRefLength*cos(nodeAngles));
  cellInfo.yPosition = bsxfun(@plus, 0, cellInfo.internalRefLength*sin(nodeAngles));
  
  cellInfo = setupInteriorNodes(cellInfo);
  cellInfo = nodeInfo(cellInfo, s);
end

%{
Barber's function, (somewhat) adapted to this program.
Function to generate the mesh/set up the interior nodes.
NOTE: Ideally, more variable names will be made more descriptive.
%}
function cellInfo = setupInteriorNodes(cellInfo)
  boungeominfo = [2; cellInfo.externalNodeCount; cellInfo.xPosition; cellInfo.yPosition];
  nameofbound = ['cell']';
  %  One can use set formulas to intersect and union bounded regions,
  %  but we just have one bounded region here
  setformulas = 'cell';
  [dl,bt] = decsg(boungeominfo,setformulas,nameofbound);
  pdegplot(dl,'EdgeLabels','on','FaceLabels','on');
  model = createpde;
  geometryFromEdges(model,dl);
  generateMesh(model,'Hmax',cellInfo.externalRefLength','Hmin',cellInfo.externalRefLength,'Hgrad',1,...
    'GeometricOrder','linear');
  pdemesh(model);
    
  tol = 1e-12;
  [a,b] = ismembertol([cellInfo.xPosition,cellInfo.yPosition],model.Mesh.Nodes',tol,'ByRows',true);
  if ~isequal(b,[1:cellInfo.externalNodeCount]')
    keyboard;
  end
    
  %  Redefine c.xn to correspond to both exterior and interior nodes
  cellInfo.xPosition = model.Mesh.Nodes(1,:)';
  cellInfo.yPosition = model.Mesh.Nodes(2,:)';
    
  %  Element matrix, used to extract reasonable interior elements...shorter
  %  name, easier to use.
  elemat = model.Mesh.Elements;

  cellInfo.totalNodeCount = numel(cellInfo.xPosition);
  cellInfo.internalNodeCount = cellInfo.totalNodeCount-cellInfo.externalNodeCount;

  %  For the "connection cells", for each connected point the first column 
  %  is its indices in terms of xs and ys, the second column is whether
  %  it is an external (false) or internal (true) node, the third column 
  %  is the index in terms of xm-ym or xn-yn depending on whichever is
  %  appropriate according to the second column.
  for nc = 1:numel(cellInfo.xPosition) % loop through total number of elements (iterating with 'nc')
    tmp = elemat(:,find(any(nc == elemat))); % tmp contains the triangles that contain an 'nc'
    tmp2 = setdiff(tmp(:),nc); % tmp2 contains a list of nodes that 'nc' is adjacent to
    if nc <= cellInfo.externalNodeCount % for all 'nc' that correspond to an external node, set 'nnext' to be the next external node 
      nnext = mod(nc,cellInfo.externalNodeCount)+1; % if we are at the last external node, set 'nnext' to the first external node
      start_ind = nnext;
    else
      start_ind = tmp2(1); % for all 'nc' that correspond to an internal node, set 'start_ind' to the first in the list of adjacent nodes
    end
    
    inds = start_ind;
    for connc = 1:numel(tmp2)-1
      tri_ind = find(any(inds(connc) == tmp),1);
      inds(connc+1) = setdiff(tmp(:,tri_ind),[nc,inds(connc)]);
      [m,n] = size(tmp);
      tmp = tmp(:,setdiff(1:n,tri_ind));
    end
    
    % makes sure numbering/orientation is correct (?)
    if ispolycw(cellInfo.xPosition(inds),cellInfo.yPosition(inds)) 
      inds = circshift(fliplr(inds),1);
    end
    
    cellInfo.nodesAdjacent{nc} = inds;
  end
  
  %  For plotting later, a list of all the indices corresponding to the
  %  line segments
  inds = [1:cellInfo.externalNodeCount]';
  cellInfo.externalLineSegments = [inds,mod(inds,cellInfo.externalNodeCount)+1];
  tmp = sort(cellInfo.externalLineSegments,2);
  cellInfo.lineSegments = [];
  for nc = 1:cellInfo.totalNodeCount
    for cnc = 1:numel(cellInfo.nodesAdjacent{nc})
      cellInfo.lineSegments = [cellInfo.lineSegments;nc,cellInfo.nodesAdjacent{nc}(cnc)];
    end
  end
  %  There are twice as many indices as we actually need for our plots,
  %  sort and eliminate duplicates:
  cellInfo.lineSegments = unique(sort(cellInfo.lineSegments,2),'rows');
  cellInfo.internalLineSegments = setdiff(cellInfo.lineSegments,tmp,'rows');

  %  Check the connectivities
  do_connectivity_plot = false;
  if do_connectivity_plot
    myc = 'rgbmyck'; myc = [myc,myc]; myc = [myc,myc]; myc = [myc,myc];
    close(figure(1)); figure(1); plot(cellInfo.xPosition,cellInfo.yPosition,'k'); hold on;
    for nc = 1:numel(cellInfo.nodesAdjacent)
      for connc = 1:numel(cellInfo.nodesAdjacent{nc})
        plot([cellInfo.xPosition(nc),cellInfo.xPosition(cellInfo.nodesAdjacent{nc}(connc))],...
          [cellInfo.yPosition(nc),cellInfo.yPosition(cellInfo.nodesAdjacent{nc}(connc))],myc(nc),'LineWidth',2);
        pause
      end
    end
  end
            
  warning(['Mesh not yet perfected (e.g. equidistant nodes, ',...
    'equal connectivity, others?']);
end

%{
Barber's function, still needs some adaptation...
Function to calculate more info about the cell
NOTE: find which fields are unnecessary
NOTE: triareainfo function?
NOTE: more descriptive variable names
%}
function cellInfo = nodeInfo(cellInfo,s)
  %  STUB
  %  Because this info is stored for each node, some info is redundant. For
  %  instance, each length will get stored twice because each edge has two
  %  nodes.  This is for simplicity and should have negligible impact on
  %  efficiency since the solving process should dominate cpu time.
  cellInfo.lengths = cellInfo.nodesAdjacent;
  cellInfo.dxs = cellInfo.nodesAdjacent;
  cellInfo.dys = cellInfo.nodesAdjacent;
  cellInfo.nxs = cellInfo.nodesAdjacent;
  cellInfo.nys = cellInfo.nodesAdjacent;
  if isfield(cellInfo,'eareas'), firsttime = false; else, firsttime = true; end
  %  These should eventually be estimated in a better way (in particular,
  %  we want a circular shape to correspond to a steady configuration) but
  %  for now we just make them 1s.
  for nc = 1:numel(cellInfo.nodesAdjacent)
    cnis = cellInfo.nodesAdjacent{nc};
    cellInfo.dxs{nc} = cellInfo.xPosition(cnis)-cellInfo.xPosition(nc);
    cellInfo.dys{nc} = cellInfo.yPosition(cnis)-cellInfo.yPosition(nc);
    cellInfo.lengths{nc} = sqrt(cellInfo.dxs{nc}.^2+cellInfo.dys{nc}.^2);
    cellInfo.txs{nc} = cellInfo.dxs{nc}./cellInfo.lengths{nc};
    cellInfo.tys{nc} = cellInfo.dys{nc}./cellInfo.lengths{nc};
    %  These make unit vectors perpendicular to the tangent vectors so that
    %  "normal" crossed with "tangent" yields 1.  The current
    %  counterclockwise orientation of external nodes makes the normal from
    %  external node i to i+1 point outwards (tangent vector points from
    %  node i to i+1)
    cellInfo.nxs{nc} = cellInfo.tys{nc};
    cellInfo.nys{nc} = -cellInfo.txs{nc};
    
%     %  Debugging
%     close(figure(5)); figure(5)
%     quiver(c.xn(nc)*ones(size(cnis)),c.yn(nc)*ones(size(cnis)),...
%       c.txs{nc}'.*c.lengths{nc}',c.tys{nc}'.*c.lengths{nc}',0);
%     hold on
%     midxs = (c.xn(nc)*ones(size(cnis))+c.xn(cnis)')./2;
%     midys = (c.yn(nc)*ones(size(cnis))+c.yn(cnis)')./2;
%     plot(c.xn(cnis),c.yn(cnis),'x');
%     quiver(midxs,midys,c.nxs{nc}',c.nys{nc}');
%     axis equal;
      
    % Between each pair of connections, there is an angle, which we
    % calculate:
    cnisp1 = circshift(1:numel(cellInfo.txs{nc}),-1); % apply a shift of elements to easily compare this node with the next node
    cellInfo.crosses{nc} = cellInfo.txs{nc}.*cellInfo.tys{nc}(cnisp1)-...
      cellInfo.tys{nc}.*cellInfo.txs{nc}(cnisp1);
    cellInfo.dots{nc} = cellInfo.txs{nc}.*cellInfo.txs{nc}(cnisp1)+...
      cellInfo.tys{nc}.*cellInfo.tys{nc}(cnisp1);
    cellInfo.alphs{nc} = atan2(cellInfo.crosses{nc},cellInfo.dots{nc});
    cellInfo.alphs{nc} = cellInfo.alphs{nc}+2*pi*(cellInfo.alphs{nc}<0);
    sum(cellInfo.alphs{nc});
    %  Between each pair of connections, there is also an element with
    %  associated info (like area of that element)
    if nc <= cellInfo.externalNodeCount
      axs = cellInfo.xPosition(cnis(1:end-1),:); ays = cellInfo.yPosition(cnis(1:end-1),:);
      bxs = cellInfo.xPosition(cnis(2:end),:); bys = cellInfo.yPosition(cnis(2:end),:);
    else
      axs = cellInfo.xPosition(cnis,:); ays = cellInfo.yPosition(cnis,:);
      bxs = cellInfo.xPosition(circshift(cnis,-1),:); bys = cellInfo.yPosition(circshift(cnis,-1),:);
    end
    xxi = cellInfo.xPosition(nc)*ones(size(axs)); xyi = cellInfo.yPosition(nc)*ones(size(bxs));
    a = [axs,ays]; b = [bxs,bys]; xi = [xxi,xyi];
    if firsttime
      cellInfo.eareas{nc} = triangleAreaInfo(xi,a,b);
    else
      %  c.anfs-"normalized forces".  This is the "geometric portion" of
      %  the area-related forces and multiplying by ka, the area elastic
      %  force modulus, will recover the actual forces
      [cellInfo.eareas{nc},cellInfo.anfs{nc}] = ...
        triareainfo(xi,a,b,cellInfo.earearefs{nc},1);
    end
  end
  %  Store area and "volume estimate"
  cellInfo.area = polyarea(cellInfo.xPosition(1:cellInfo.externalNodeCount),cellInfo.yPosition(1:cellInfo.externalNodeCount));
  cellInfo.volest = "unimplemented for now";
  %est_3d_volume(cellInfo.xPosition(1:cellInfo.externalNodeCount),...
    %cellInfo.yPosition(1:cellInfo.externalNodeCount),s.vol_est_type);
end

