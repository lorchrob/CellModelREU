%{
Main function to initialize the system
%}
function cellInfo = initializeNetwork(externalNodeCount, simulationType)  
  if exist('simulationType', 'var') && strcmp(simulationType, 'wall') 
    cellInfo.xw =   [-190, -150,  150, 150, 190, 190, 150, 150, -150, -190];
    cellInfo.yw = 2*[ -15, -2.5, -2.5, -15, -15,  15,  15, 2.5,  2.5,   15];
  else
    cellInfo.xw = [-500,  500, 500, -500]; % defaults to a square with length 1000
    cellInfo.yw = [-500, -500, 500,  500];
  end

  % set initial values for some fields of the 'cellInfo' struct
  cellInfo.externalNodeCount = externalNodeCount;
  cellInfo.radius = 12.4/2; 
  cellInfo.refArea = pi*cellInfo.radius^2;
  % underscores represent subscripts
  cellInfo.k_te = 1200; %1200
  cellInfo.k_ti = 1200; %1200
  cellInfo.mu_e = 200; %200
  cellInfo.mu_i = 100; %100
  cellInfo.k_be = 90; %90
  cellInfo.k_bi = 0;
  
  % distance from center to the external nodes
  cellInfo.internalRefLength = sqrt(2*pi*cellInfo.radius^2/cellInfo.externalNodeCount/sin(2*pi/cellInfo.externalNodeCount)); 
  cellInfo.externalRefLength = 2*sin(2*pi/cellInfo.externalNodeCount/2)*cellInfo.internalRefLength;
  
  % set up outer circle
  nodeAngles = linspace(0, 2*pi, cellInfo.externalNodeCount+1)';
  nodeAngles = nodeAngles(1:end-1);
  cellInfo.xPosition = bsxfun(@plus, 0, cellInfo.internalRefLength*cos(nodeAngles));
  cellInfo.yPosition = bsxfun(@plus, 0, cellInfo.internalRefLength*sin(nodeAngles));
  
  cellInfo = setupInteriorNodes(cellInfo);
  cellInfo = calculateNodeInfo(cellInfo);
  
  %  For plotting later, a list of all the indices corresponding to the
  %  line segments
  inds = [1:cellInfo.externalNodeCount]';
  cellInfo.externalLineSegments = [inds,mod(inds,cellInfo.externalNodeCount)+1];
  tmp = sort(cellInfo.externalLineSegments,2);
  cellInfo.lineSegments = [];
  for i = 1:cellInfo.totalNodeCount
    for cnc = 1:numel(cellInfo.nodesAdjacent{i})
      cellInfo.lineSegments = [cellInfo.lineSegments;i,cellInfo.nodesAdjacent{i}(cnc)];
    end
  end
  %  There are twice as many indices as we actually need for our plots,
  %  sort and eliminate duplicates:
  cellInfo.lineSegments = unique(sort(cellInfo.lineSegments,2),'rows');
  cellInfo.internalLineSegments = setdiff(cellInfo.lineSegments,tmp,'rows');
  
  % For now, we say that no nodes are fixed and there are no external
  % forces (may be updated later)
  cellInfo.isFixed = false(cellInfo.totalNodeCount, 1);
  cellInfo.externalForces = zeros(cellInfo.totalNodeCount, 2);
  cellInfo.xwf = zeros(cellInfo.totalNodeCount,1);
  cellInfo.ywf = cellInfo.xwf;
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
  model = createpde;
  geometryFromEdges(model, dl);
  generateMesh(model, 'Hmax', cellInfo.externalRefLength', 'Hmin', cellInfo.externalRefLength, 'Hgrad', 1,...
    'GeometricOrder', 'linear');
    
  tol = 1e-12;
  [a,b] = ismembertol([cellInfo.xPosition, cellInfo.yPosition], model.Mesh.Nodes', tol, 'ByRows', true);
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
  for i = 1:numel(cellInfo.xPosition) % loop through total number of elements (iterating with 'i')
    tmp = elemat(:,find(any(i == elemat))); % tmp contains the triangles that contain an 'i'
    tmp2 = setdiff(tmp(:),i); % tmp2 contains a list of nodes that 'i' is adjacent to
    if i <= cellInfo.externalNodeCount % for all 'i' that correspond to an external node, set 'nnext' to be the next external node 
      nnext = mod(i, cellInfo.externalNodeCount)+1; % if we are at the last external node, set 'nnext' to the first external node
      start_ind = nnext;
    else
      start_ind = tmp2(1); % for all 'i' that correspond to an internal node, set 'start_ind' to the first in the list of adjacent nodes
    end
    
    inds = start_ind;
    for connc = 1:numel(tmp2)-1
      tri_ind = find(any(inds(connc) == tmp),1);
      inds(connc+1) = setdiff(tmp(:,tri_ind),[i,inds(connc)]);
      [m,n] = size(tmp);
      tmp = tmp(:,setdiff(1:n,tri_ind));
    end
    
    % makes sure numbering/orientation is correct (?)
    if ispolycw(cellInfo.xPosition(inds),cellInfo.yPosition(inds)) 
      inds = circshift(fliplr(inds),1);
    end
    
    cellInfo.nodesAdjacent{i} = inds;
  end
end



