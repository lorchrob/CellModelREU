%{
Function to calculate the forces on nodes due to a wall. Should be run at
every timestep.

Function takes in:
 - an x 'wall'
 - cellInfo

Function then
 - calculates the forces on each node due to the wall
 - sets these forces into 'cellInfo.externalForces'

Function returns
 - 'cellInfo' with updated external forces


xw =   [-190, -150,  150, 150, 190, 190, 150, 150, -150, -190]; %#ok<UNRCH>
yw = 2*[ -15, -2.5, -2.5, -15, -15,  15,  15, 2.5,  2.5,   15];

%}

function cellInfoNew = updateForceCellWall(cellInfo, xw, yw)
  cellInfoNew = cellInfo;
%  xw =   [-190, -150,  150, 150, 190, 190, 150, 150, -150, -190];
%  yw = 2*[ -15, -2.5, -2.5, -15, -15,  15,  15, 2.5,  2.5,   15];
  wallForceMag = 2000; % Magnitude of force from wall
  

  [dists,dot_dist,overlap,norm_x,norm_y,norm_arc,type] = ...
    dist_from_pt_to_line_segs(cellInfo.xPosition, cellInfo.yPosition, xw, yw);
  %  This selects out only the segments on the boundary that the cell nodes
  %  are closest to: mindists-vector of n_nodes with distances and
  %  inds-vector of indices corresponding to the line segment on the
  %  boundary that is closest to each node
  [mindists,inds] = min(dists,[],1);
  %  Convert these indices into "linear indices" for easier use with
  %  matrices
  lininds = sub2ind(size(dists),inds,1:cellInfo.totalNodeCount);
  %  Decide if the cell nodes are inside (1) or outside (0) of the boundary
  my_eps = 0.1;
  in_or_out = inpolygon(cellInfo.xPosition,cellInfo.yPosition,xw,yw);

  ramp_func = @(x,e) (x > -e).*(x < 0).*(x+e).^2.*(e-2*x)./e.^3+(x >= 0);
  
  % Wall forces
  fxnwf = (wallForceMag*(norm_x(lininds).*...
    ramp_func((1-2*in_or_out').*mindists,my_eps)))';
  fynwf = (wallForceMag*(norm_y(lininds).*...
    ramp_func((1-2*in_or_out').*mindists,my_eps)))';
  
  cellInfoNew.xwf(:) = fxnwf;
  cellInfoNew.ywf(:) = fynwf;



%   for i = 1:cellInfo.totalNodeCount
%     if cellInfoNew.yPosition(i) - max(xw) > 0
%       cellInfoNew.externalForces(i,2) = -abs(cellInfoNew.yPosition(i) - max(xw)) * 200;
%     end
%     
%     if min(xw) - cellInfoNew.yPosition(i) > 0
%       cellInfoNew.externalForces(i,2) = abs(min(xw) - cellInfoNew.yPosition(i)) * 200;
%     end
%   end
  
end


