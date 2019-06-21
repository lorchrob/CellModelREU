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
  cellInfo.externalNodeCount = externalNodeCount;
  cellInfo.refLength = 3; % placeholder
  cellInfo.refArea = 3; % placeholder
  cellInfo.tensions = []; % ? maybe ?
  cellInfo.lengths = []; 
  
  % used ONLY for determining internalRefLength (from Barber's code)
  cellInfo.lengthScale = 12.4/2; 
  
  % distance from center to the external nodes
  cellInfo.internalRefLength = sqrt(2*pi*cellInfo.lengthScale^2/cellInfo.externalNodeCount/sin(2*pi/cellInfo.externalNodeCount)); 
  
  % set up outer circle
  nodeAngles = linspace(0, 2*pi, cellInfo.externalNodeCount+1)';
  nodeAngles = nodeAngles(1:end-1);
  cellInfo.xPosition = bsxfun(@plus, 0, cellInfo.internalRefLength*cos(nodeAngles));
  cellInfo.yPosition = bsxfun(@plus, 0, cellInfo.internalRefLength*sin(nodeAngles));
  
  % complete initialization of the 'cellInfo' struct
  % NOTE: still unsure of return type of 'createMesh()' function
  cellInfo = createMesh(cellInfo);
  plot(cellInfo.xPosition, cellInfo.yPosition);
end

%{
Function to generate the triangular mesh for the cell interior
NOTE: still unsure of return type for this function
%}
function cellInfo = createMesh(cellInfo)
  % STUB
end

%{
Barber's function, needs some adaptation...
%}
function cellInfo = setup_int_nodes(cellInfo)
  % STUB
  if cellInfo.mult_int_nodes
    boungeominfo = [2; cellInfo.n_enodes; cellInfo.xn; cellInfo.yn];
    nameofbound = ['cell']';
    %  One can use set formulas to intersect and union bounded regions,
    %  but we just have one bounded region here
    setformulas = 'cell';
    [dl,bt] = decsg(boungeominfo,setformulas,nameofbound);
    pdegplot(dl,'EdgeLabels','on','FaceLabels','on');
    model = createpde;
    geometryFromEdges(model,dl);
    generateMesh(model,'Hmax',cellInfo.loref','Hmin',cellInfo.loref,'Hgrad',1,...
      'GeometricOrder','linear');
    pdemesh(model);
    
    tol = 1e-12;
    [a,b] = ismembertol([cellInfo.xn,cellInfo.yn],model.Mesh.Nodes',tol,'ByRows',true);
    if ~isequal(b,[1:cellInfo.n_enodes]')
      keyboard;
    end
    
    %  Redefine c.xn to correspond to both exterior and interior nodes
    cellInfo.xn = model.Mesh.Nodes(1,:)';
    cellInfo.yn = model.Mesh.Nodes(2,:)';
    
    %  Element matrix, used to extract reasonable interior elements...shorter
    %  name, easier to use.
    elemat = model.Mesh.Elements;
  else
    cellInfo.xn = [cellInfo.xn;cellInfo.xm];
    cellInfo.yn = [cellInfo.yn;cellInfo.ym];
    
    for enc = 1:(cellInfo.n_enodes-1)
      elemat(:,enc) = [enc;cellInfo.n_enodes+1;enc+1];
    end
    elemat(:,cellInfo.n_enodes) = [cellInfo.n_enodes;cellInfo.n_enodes+1;1];
  end
  cellInfo.n_nodes = numel(cellInfo.xn);
  cellInfo.n_inodes = cellInfo.n_nodes-cellInfo.n_enodes;

  %  For the "connection cells", for each connected point the first column 
  %  is its indices in terms of xs and ys, the second column is whether
  %  it is an external (false) or internal (true) node, the third column 
  %  is the index in terms of xm-ym or xn-yn depending on whichever is
  %  appropriate according to the second column.
  for nc = 1:numel(cellInfo.xn)
    tmp = elemat(:,find(any(nc == elemat)));
    tmp2 = setdiff(tmp(:),nc);
    if nc <= cellInfo.n_enodes
      nnext = mod(nc,cellInfo.n_enodes)+1;
      start_ind = nnext;
    else
      start_ind = tmp2(1);
    end
    inds = start_ind;
    for connc = 1:numel(tmp2)-1
      tri_ind = find(any(inds(connc) == tmp),1);
      inds(connc+1) = setdiff(tmp(:,tri_ind),[nc,inds(connc)]);
      [m,n] = size(tmp);
      tmp = tmp(:,setdiff(1:n,tri_ind));
    end
    if ispolycw(cellInfo.xn(inds),cellInfo.yn(inds))
      inds = circshift(fliplr(inds),1);
    end
    cellInfo.nc{nc} = inds;
  end
  
  %  For plotting later, a list of all the indices corresponding to the
  %  line segments
  inds = [1:cellInfo.n_enodes]';
  cellInfo.enlist = [inds,mod(inds,cellInfo.n_enodes)+1];
  tmp = sort(cellInfo.enlist,2);
  cellInfo.nlist = [];
  for nc = 1:cellInfo.n_nodes
    for cnc = 1:numel(cellInfo.nc{nc})
      cellInfo.nlist = [cellInfo.nlist;nc,cellInfo.nc{nc}(cnc)];
    end
  end
  %  There are twice as many indices as we actually need for our plots,
  %  sort and eliminate duplicates:
  cellInfo.nlist = unique(sort(cellInfo.nlist,2),'rows');
  cellInfo.inlist = setdiff(cellInfo.nlist,tmp,'rows');

  %  Check the connectivities
  do_connectivity_plot = false;
  if do_connectivity_plot
    myc = 'rgbmyck'; myc = [myc,myc]; myc = [myc,myc]; myc = [myc,myc];
    close(figure(1)); figure(1); plot(cellInfo.xn,cellInfo.yn,'k'); hold on;
    for nc = 1:numel(cellInfo.nc)
      for connc = 1:numel(cellInfo.nc{nc})
        plot([cellInfo.xn(nc),cellInfo.xn(cellInfo.nc{nc}(connc))],...
          [cellInfo.yn(nc),cellInfo.yn(cellInfo.nc{nc}(connc))],myc(nc),'LineWidth',2);
        pause
      end
    end
  end
            
  warning(['Mesh not yet perfected (e.g. equidistant nodes, ',...
    'equal connectivity, others?']);
end

%{
Barber's function, still needs some adaptation...
%}
function cellInfo = node_info(cellInfo,s)
  %  STUB
  %  Because this info is stored for each node, some info is redundant. For
  %  instance, each length will get stored twice because each edge has two
  %  nodes.  This is for simplicity and should have negligible impact on
  %  efficiency since the solving process should dominate cpu time.
  cellInfo.ls = cellInfo.nc;
  cellInfo.dxs = cellInfo.nc;
  cellInfo.dys = cellInfo.nc;
  cellInfo.nxs = cellInfo.nc;
  cellInfo.nys = cellInfo.nc;
  if isfield(cellInfo,'eareas'), firsttime = false; else, firsttime = true; end
  %  These should eventually be estimated in a better way (in particular,
  %  we want a circular shape to correspond to a steady configuration) but
  %  for now we just make them 1s.
  for nc = 1:numel(cellInfo.nc)
    cnis = cellInfo.nc{nc};
    cellInfo.dxs{nc} = cellInfo.xn(cnis)-cellInfo.xn(nc);
    cellInfo.dys{nc} = cellInfo.yn(cnis)-cellInfo.yn(nc);
    cellInfo.ls{nc} = sqrt(cellInfo.dxs{nc}.^2+cellInfo.dys{nc}.^2);
    cellInfo.txs{nc} = cellInfo.dxs{nc}./cellInfo.ls{nc};
    cellInfo.tys{nc} = cellInfo.dys{nc}./cellInfo.ls{nc};
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
%       c.txs{nc}'.*c.ls{nc}',c.tys{nc}'.*c.ls{nc}',0);
%     hold on
%     midxs = (c.xn(nc)*ones(size(cnis))+c.xn(cnis)')./2;
%     midys = (c.yn(nc)*ones(size(cnis))+c.yn(cnis)')./2;
%     plot(c.xn(cnis),c.yn(cnis),'x');
%     quiver(midxs,midys,c.nxs{nc}',c.nys{nc}');
%     axis equal;
      
    % Between each pair of connections, there is an angle, which we
    % calculate:
    cnisp1 = circshift(1:numel(cellInfo.txs{nc}),-1);
    cellInfo.crosses{nc} = cellInfo.txs{nc}.*cellInfo.tys{nc}(cnisp1)-...
      cellInfo.tys{nc}.*cellInfo.txs{nc}(cnisp1);
    cellInfo.dots{nc} = cellInfo.txs{nc}.*cellInfo.txs{nc}(cnisp1)+...
      cellInfo.tys{nc}.*cellInfo.tys{nc}(cnisp1);
    cellInfo.alphs{nc} = atan2(cellInfo.crosses{nc},cellInfo.dots{nc});
    cellInfo.alphs{nc} = cellInfo.alphs{nc}+2*pi*(cellInfo.alphs{nc}<0);
    sum(cellInfo.alphs{nc});
    %  Between each pair of connections, there is also an element with
    %  associated info (like area of that element)
    if nc <= cellInfo.n_enodes
      axs = cellInfo.xn(cnis(1:end-1),:); ays = cellInfo.yn(cnis(1:end-1),:);
      bxs = cellInfo.xn(cnis(2:end),:); bys = cellInfo.yn(cnis(2:end),:);
    else
      axs = cellInfo.xn(cnis,:); ays = cellInfo.yn(cnis,:);
      bxs = cellInfo.xn(circshift(cnis,-1),:); bys = cellInfo.yn(circshift(cnis,-1),:);
    end
    xxi = cellInfo.xn(nc)*ones(size(axs)); xyi = cellInfo.yn(nc)*ones(size(bxs));
    a = [axs,ays]; b = [bxs,bys]; xi = [xxi,xyi];
    if firsttime
      cellInfo.eareas{nc} = triareainfo(xi,a,b);
    else
      %  c.anfs-"normalized forces".  This is the "geometric portion" of
      %  the area-related forces and multiplying by ka, the area elastic
      %  force modulus, will recover the actual forces
      [cellInfo.eareas{nc},cellInfo.anfs{nc}] = ...
        triareainfo(xi,a,b,cellInfo.earearefs{nc},1);
    end
  end
  %  Store area and "volume estimate"
  cellInfo.area = polyarea(cellInfo.xn(1:cellInfo.n_enodes),cellInfo.yn(1:cellInfo.n_enodes));
  cellInfo.volest = est_3d_volume(cellInfo.xn(1:cellInfo.n_enodes),...
    cellInfo.yn(1:cellInfo.n_enodes),s.vol_est_type);
end

