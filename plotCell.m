%{
Function to plot a cell with a red membrane and a blue interior
%}
function plotCell(cellInfo)

  minX = min(cellInfo.xPosition);
  maxX = max(cellInfo.xPosition);
  maxY = max(cellInfo.yPosition);
  minY = min(cellInfo.yPosition);
  plot([cellInfo.xPosition(cellInfo.internalLineSegments(1:end, 1)), cellInfo.xPosition(cellInfo.internalLineSegments(1:end, 2))]',... 
       [cellInfo.yPosition(cellInfo.internalLineSegments(1:end, 1)), cellInfo.yPosition(cellInfo.internalLineSegments(1:end, 2))]',...
       'b');
  a = 15;
%   axis([-a a -a a])
  axis([minX-5, maxX + 5, minY-5, maxY+5])
  hold on
  plot([cellInfo.xPosition(cellInfo.externalLineSegments(1:end, 1)), cellInfo.xPosition(cellInfo.externalLineSegments(1:end, 2))]',... 
       [cellInfo.yPosition(cellInfo.externalLineSegments(1:end, 1)), cellInfo.yPosition(cellInfo.externalLineSegments(1:end, 2))]',...
       'r');
  hold on
   
  plot([ cellInfo.xw(1:end-1) ; cellInfo.xw(2:end) ]',...
       [ cellInfo.yw(1:end-1); cellInfo.yw(2:end) ]',...
       'm');
  axis manual
  
%   if isfield(cellInfo, "modelType") && cellInfo.modelType == "timeStepper"
%     quiver(cellInfo.xPosition(:), cellInfo.yPosition(:),...
%            cellInfo.xVelocity(:), cellInfo.yVelocity(:),...
%            0, 'k');
%   end
  
  hold off
end