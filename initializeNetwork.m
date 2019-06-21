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
NOTE: not sure yet if necessary to use all three returned structs
%}
function [cellInfo, systemInfo, programInfo] = initializeNetwork(nodeCount)  
  % set initial values for all fields of the 'cellInfo' struct
  % NOTE: still unsure which values are given and which are calculated
  cellInfo.nodeCount = nodeCount;
  cellInfo.referenceLength = 3; % placeholder
  cellInfo.referenceArea = 3; % placeholder
  cellInfo.tensions = []; % initialization, will contain all the tension forces
  cellInfo.distances = []; % initialization, will contain the distance of every connection
  
  % complete initialization of the 'cellInfo' struct
  % NOTE: still unsure of return type of 'createMesh()' function
  cellInfo = createMesh(cellInfo);
  
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
function c = setup_int_nodes(c)

  if c.mult_int_nodes
    boungeominfo = [2;c.n_enodes;c.xn;c.yn];
    nameofbound = ['cell']';
    %  One can use set formulas to intersect and union bounded regions,
    %  but we just have one bounded region here
    setformulas = 'cell';
    [dl,bt] = decsg(boungeominfo,setformulas,nameofbound);
    pdegplot(dl,'EdgeLabels','on','FaceLabels','on');
    model = createpde;
    geometryFromEdges(model,dl);
    generateMesh(model,'Hmax',c.loref','Hmin',c.loref,'Hgrad',1,...
      'GeometricOrder','linear');
    pdemesh(model);
    
    tol = 1e-12;
    [a,b] = ismembertol([c.xn,c.yn],model.Mesh.Nodes',tol,'ByRows',true);
    if ~isequal(b,[1:c.n_enodes]')
      keyboard;
    end
    
    %  Redefine c.xn to correspond to both exterior and interior nodes
    c.xn = model.Mesh.Nodes(1,:)';
    c.yn = model.Mesh.Nodes(2,:)';
    
    %  Element matrix, used to extract reasonable interior elements...shorter
    %  name, easier to use.
    elemat = model.Mesh.Elements;
  else
    c.xn = [c.xn;c.xm];
    c.yn = [c.yn;c.ym];
    
    for enc = 1:(c.n_enodes-1)
      elemat(:,enc) = [enc;c.n_enodes+1;enc+1];
    end
    elemat(:,c.n_enodes) = [c.n_enodes;c.n_enodes+1;1];
  end
  c.n_nodes = numel(c.xn);
  c.n_inodes = c.n_nodes-c.n_enodes;

  %  For the "connection cells", for each connected point the first column 
  %  is its indices in terms of xs and ys, the second column is whether
  %  it is an external (false) or internal (true) node, the third column 
  %  is the index in terms of xm-ym or xn-yn depending on whichever is
  %  appropriate according to the second column.
  for nc = 1:numel(c.xn)
    tmp = elemat(:,find(any(nc == elemat)));
    tmp2 = setdiff(tmp(:),nc);
    if nc <= c.n_enodes
      nnext = mod(nc,c.n_enodes)+1;
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
    if ispolycw(c.xn(inds),c.yn(inds))
      inds = circshift(fliplr(inds),1);
    end
    c.nc{nc} = inds;
  end
  
  %  For plotting later, a list of all the indices corresponding to the
  %  line segments
  inds = [1:c.n_enodes]';
  c.enlist = [inds,mod(inds,c.n_enodes)+1];
  tmp = sort(c.enlist,2);
  c.nlist = [];
  for nc = 1:c.n_nodes
    for cnc = 1:numel(c.nc{nc})
      c.nlist = [c.nlist;nc,c.nc{nc}(cnc)];
    end
  end
  %  There are twice as many indices as we actually need for our plots,
  %  sort and eliminate duplicates:
  c.nlist = unique(sort(c.nlist,2),'rows');
  c.inlist = setdiff(c.nlist,tmp,'rows');

  %  Check the connectivities
  do_connectivity_plot = false;
  if do_connectivity_plot
    myc = 'rgbmyck'; myc = [myc,myc]; myc = [myc,myc]; myc = [myc,myc];
    close(figure(1)); figure(1); plot(c.xn,c.yn,'k'); hold on;
    for nc = 1:numel(c.nc)
      for connc = 1:numel(c.nc{nc})
        plot([c.xn(nc),c.xn(c.nc{nc}(connc))],...
          [c.yn(nc),c.yn(c.nc{nc}(connc))],myc(nc),'LineWidth',2);
        pause
      end
    end
  end
            
  warning(['Mesh not yet perfected (e.g. equidistant nodes, ',...
    'equal connectivity, others?']);

end

