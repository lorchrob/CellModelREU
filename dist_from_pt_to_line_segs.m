%{
Function from Barber. Helper function to help calculate wall force.
%}
function [dists,dot_dist,overlap,norm_x,norm_y,norm_arc,type] = ...
  dist_from_pt_to_line_segs(pts_xs,pts_ys,lseg_xs,lseg_ys)

  %  Finds minimum distances from each point in a list of points to a set
  %    of line segments including additional information
  %  The primary algorithm is for each point-line segment pairing
  %    Find the distance from the point to both endpoints of the line 
  %    Find the distance from the point to the line that passes through the
  %      line segment (using the dot product)
  %    Use dot product information to select which of the three distances
  %      actually corresponds to the minimum point-line segment distance
  %
  %  Inputs:
  %  pts_xs,pts_ys-x and y coordinates for the list of points
  %  lseg_xs,lseg_ys-the endpoints for the line segments in question.  We
  %    assume that this a connected set of line segments where lseg_?s(1)
  %    to lseg_?s(2) is the first line segment, lseg_?s(2) to lseg_?s(3) is
  %    the second line segment, and so on until the end of the list.
  %
  %  Outputs:
  %  In all cases the outputs below are matrices of size (numel(lseg_xs)-1)
  %    by numel(pts_xs) for a matrix M, M(i,j) corresponds to the ith
  %    point and its relationship to the jth line segment.  Call the ith
  %    point P and the endpoints of the jth line segment Q and R.  Also
  %    call the point on QR that lies closest to P point S.
  %  dists-minimum distance from the point to the line segment (|PS|)
  %  dot_dist-scalar projection of QP onto QR
  %  overlap-Minimum of |QS| and |SR|.  Note that in many cases overlap =
  %    0.  This is used for calculating lubrication forces later on.
  %  norm_x/norm_y-unit vector point from S to P.  There is a sign
  %    ambiguity here that shouldn't matter for our earliest models.
  %  norm_arc-when S corresponds to Q or R, this is just 0 or 1.  
  %    Otherwise, this is dot_dist/length of line segment, i.e. a
  %    "normalized arclength".
  %  type: 0-PS is perpendicular to QR
  %    -1-PS not perp to QR and P is closest to left endpoint
  %    1-PS not perp to QR and P i closest to right endpoint
  
  %  Initial setup
  %  Number of points in the point vector
  n_pts = numel(pts_xs);
  %  Number of nodes in the line segment vector...note that this
  %  means that there are n_nod-1 line segments
  n_lseg = numel(lseg_xs);  
  %  To make best use of bsxfun, we reshape points to correspond to row
  %  vectors and line segment points to correspond to column vectors
  pts_xs = reshape(pts_xs,1,n_pts);
  pts_ys = reshape(pts_ys,1,n_pts);
  lseg_xs = reshape(lseg_xs,n_lseg,1);
  lseg_ys = reshape(lseg_ys,n_lseg,1);  

  %  Find minimum distances from points to line segment endpoints 
  %  endpt_d?(i,j)-n_nodxn_pts, the x and y components of the vectors from
  %  the ith line segment node to the jth point.
  endpt_dx = bsxfun(@minus,pts_xs,lseg_xs);
  endpt_dy = bsxfun(@minus,pts_ys,lseg_ys);
  %  endpt_dists(i,j)-n_nodxn_pts, the distance from node i to point j
  endpt_dist = sqrt(endpt_dx.^2+endpt_dy.^2);
  %  dist_?????(i,j)-(n_nod-1)xn_pts, the distance from the lower and upper
  %  endpoints for the ith line segment to point j
  dist_lower = endpt_dist(1:end-1,:);
  dist_upper = endpt_dist(2:end,:);
  %  norm_endpt_d?(i,j)-n_ptsxn_nod, the unit-vector that points from the
  %  ith line segment node to the jth point
  %  tan_endpt_d?(i,j)-n_ptsxn_nod, a unit-vector that is perpendicular to
  %  the above normendptd?
  %  Corresponding "normal" and "tangent" vectors for special use only when
  %  two corners are in close proximity (as opposed to the usual point
  %  being close to a line segment).  We note that there is an ambiguity in
  %  the definition of the tangent but that is ok because the sign doesn't
  %  matter as the expression for the force (2.19 in my final dissertation)
  %  does not depend on the signs of the normal and tangent vectors.
  norm_endpt_dx = endpt_dx./endpt_dist;
  norm_endpt_dy = endpt_dy./endpt_dist;
  tan_endpt_dx = norm_endpt_dy;
  tan_endpt_dy = -norm_endpt_dx;
  
  %  Find distance from the points to the line that passes through each
  %  line segment
  %  lin_d?(j)-(n_nod-1) The vector corresponding to each line segment
  lin_dx = diff(lseg_xs);
  lin_dy = diff(lseg_ys);
  %  lin_dist(j)-(n_nod-1) Length of the vector corresponding to each line
  %  segment.
  lin_dist = sqrt(lin_dx.^2+lin_dy.^2);
  %  tan_d?(j)-(n_nod-1) The tangent vectors (unit length) along each line
  %  segment
  tan_dx = lin_dx./lin_dist;
  tan_dy = lin_dy./lin_dist;
  %  norm_d?_o(j)-(n_nod-1), corresponding normal vectors for above tangent
  %    vectors...as already mentioned we don't worry about the sign
  %    ambiguity
  norm_dx_o = -tan_dy;
  norm_dy_o = tan_dx;
  %  endpt_1_to_pt_?(i,j)-n_ptsx(n_nod-1) Vector from the jth node/first line
  %  segment endpoint for the jth line segment to the ith point
  endpt_1_to_pt_x = bsxfun(@minus,pts_xs,lseg_xs(1:end-1));
  endpt_1_to_pt_y = bsxfun(@minus,pts_ys,lseg_ys(1:end-1));
  %  dot_dist(i,j)-n_ptsx(n_nod-1) signed length of the projection of the
  %  ith endpt_1_to_pt vector onto the jth line segment
  dot_dist = bsxfun(@times,tan_dx,endpt_1_to_pt_x)+...
    bsxfun(@times,tan_dy,endpt_1_to_pt_y);
  %  ortho_dists(i,j)-n_ptsx(n_nod-1) signed distance from the point to the
  %  line that passes through the line segment in question
  ortho_dists = bsxfun(@times,tan_dx,endpt_1_to_pt_y)-...
    bsxfun(@times,tan_dy,endpt_1_to_pt_x);
  
  %  Decide which of the two endpoint distances is shortest and store
  %  corresponding "normal vector" information
  dists = min(dist_lower,dist_upper);
  type = (dists == dist_upper)-(dists == dist_lower);
  norm_x = (dists == dist_lower).*norm_endpt_dx(1:end-1,:)+...
    (dists == dist_upper).*norm_endpt_dx(2:end,:);
  norm_y = (dists == dist_lower).*norm_endpt_dy(1:end-1,:)+...
    (dists == dist_upper).*norm_endpt_dy(2:end,:);
  
  %  Decide whether or not to use the orthogonal distance
  %  How much the vector from first endpoint to the point overlaps the
  %  vector corresponding to the line segment.
  overlap = max(0,min(dot_dist,bsxfun(@minus,lin_dist,dot_dist)));
  %  This "normalizes" the projection distance by dividing by the total
  %  distance of the line segment.  "Normalized arclength".
  norm_arc = min(1,max(0,bsxfun(@rdivide,dot_dist,lin_dist)));
  %  Indices where the point lies closest the line segment between the two
  %  endpoints rather than to just one or the other endpoint.
  inds = (norm_arc > 0) & (norm_arc < 1) & (dists ~= 0);
  %  Correct the minimum distance estimate in those regions
  dists(inds) = abs(ortho_dists(inds));
  %  Correct the corresponding normal vectors as well.
  tmpnormdxo = norm_dx_o*ones(1,n_pts);
  tmpnormdyo = norm_dy_o*ones(1,n_pts);
  norm_x(inds) = tmpnormdxo(inds);
  norm_y(inds) = tmpnormdyo(inds);
  type(inds) = 0;
  
end