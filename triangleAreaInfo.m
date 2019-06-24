%{
This .m is copied from Barber
%}
function [area,F] = triareainfo(xi,a,b,prevareas,ka)
  %  Returns various information for an elastic triangle
  %  The "areal elastic energy" of the triangle per unit area is given as a
  %  function of "areal strain":  Ea/Ae0 = ka/2*((Ae-Ae0)/Ae0)^2.
  %  
  %  xi, a, and b are all n X (2 or 3) matrices where n is the number of
  %    triangular elements we plan to loop through
  %  xi-vertex whose force we are interested in
  %  a and b-other two vertices for a particular element
  %  previnfo-contains reference area info for use in calculating elastic
  %    forces--if empty no elastic forces are calculated, just areas
  %  ka-elastic modulus of sorts
  
%   nd = size(xi,2);
%   if nd == 2
%   %  For debugging (2d only...polyarea doesn't work for 3d polygons)
%   xs = [xi(:,1),a(:,1),b(:,1)]'; ys = [xi(:,2),a(:,2),b(:,2)]';
% %   mycols = 'rgbmcyk'; mycols = [mycols,mycols];
% %   p1 = patch(xs,ys,(1:size(xi,1))');
%   area1 = polyarea(xs,ys);
%   for j = 1:size(xi,1)
%     area2(j) = polyarea([xi(j,1),a(j,1),b(j,1)],[xi(j,2),a(j,2),b(j,2)]);
%   end
%     xi(:,3) = 0; a(:,3) = 0; b(:,3) = 0;
%   end
  v1 = a-xi; v2 = b-xi;
%   area3 = cross(v1,v2,2)./2; area3 = sqrt(sum(area3.^2,2));
%   %  Alternative area using "Lagrange's identity"
%   area4 = sqrt((sum(v1.^2,2).*sum(v2.^2,2)-sum(v1.*v2,2).^2)./4);
%   fprintf('area1 = '); fprintf('%g ',area1); fprintf('\n');
%   fprintf('area2 = '); fprintf('%g ',area2); fprintf('\n');
%   fprintf('area3 = '); fprintf('%g ',area3); fprintf('\n');
%   fprintf('area4 = '); fprintf('%g ',area4); fprintf('\n');
%   area = area4;

  area = sqrt((sum(v1.^2,2).*sum(v2.^2,2)-sum(v1.*v2,2).^2)./4);
  
  %  If there is no previous info, then we don't have a reference area and
  %  can't really calculate an elastic force
  if nargin > 3
    F = (-ka/4./area).*(area./prevareas-1).*...
      (-sum(v2.*(v2-v1),2).*v1+sum(v1.*(v2-v1),2).*v2);
    
%     dA2dl = 1/2*(-v1.*sum(v2.*v2,2)-v2.*sum(v1.*v1,2)+...
%       (v1+v2).*sum(v1.*v2,2))
%     F2 = ka*(area./prevareas-1).*(dA2dl./2./area)
%     my_bump = sqrt(eps(norm(xi)));
%     %  Debugging purposes only (2d test)
%     for elec = 1:size(xi,1)
%       area_baseline = area_dot_formula(xi(elec,:),a(elec,:),b(elec,:))^2;
%       area_bplusdx = area_dot_formula(xi(elec,:)+my_bump*[1,0,0],...
%         a(elec,:),b(elec,:))^2;
%       area_bplusdy = area_dot_formula(xi(elec,:)+my_bump*[0,1,0],...
%         a(elec,:),b(elec,:))^2;
%       darea2dl(elec,1) = (area_bplusdx-area_baseline)./my_bump;
%       darea2dl(elec,2) = (area_bplusdy-area_baseline)./my_bump;
%       energy_baseline = ...
%         energy(xi(elec,:),a(elec,:),b(elec,:),prevareas(elec),ka);
%       energy_bplusdx = ...
%         energy(xi(elec,:)+my_bump*[1,0,0],a(elec,:),b(elec,:),...
%         prevareas(elec),ka);
%       energy_bplusdy = ...
%         energy(xi(elec,:)+my_bump*[0,1,0],a(elec,:),b(elec,:),...
%         prevareas(elec),ka);
%       denergydl(elec,1) = (energy_bplusdx-energy_baseline)./my_bump;
%       denergydl(elec,2) = (energy_bplusdy-energy_baseline)./my_bump;
%     end
%     keyboard
  else
    F = 0;
  end
  
end

function area = area_dot_formula(xi,a,b)
  v1 = a-xi; v2 = b-xi;
  area = sqrt((sum(v1.^2,2).*sum(v2.^2,2)-sum(v1.*v2,2).^2)./4);
end

function energy = energy(xi,a,b,Ae0,ka)
  energy = (ka*Ae0/2)*(area_dot_formula(xi,a,b)-Ae0)^2/Ae0^2;
end