%{
Function to plot a cell. To change colors, change the last parameter to the
'plot' calls.
 %}
function plotCell(cellInfo)

  minX = min(cellInfo.xPosition);
  maxX = max(cellInfo.xPosition);
  maxY = max(cellInfo.yPosition);
  minY = min(cellInfo.yPosition);
  plot([cellInfo.xPosition(cellInfo.internalLineSegments(1:end, 1)), cellInfo.xPosition(cellInfo.internalLineSegments(1:end, 2))]',... 
       [cellInfo.yPosition(cellInfo.internalLineSegments(1:end, 1)), cellInfo.yPosition(cellInfo.internalLineSegments(1:end, 2))]',...
       'r');
  a = 15;
  %axis([-a a -a a])
  axis([minX-5, maxX + 5, minY-5, maxY+5])
  %axis([-20, 20, -8, 8]);
  
  hold on
  plot([cellInfo.xPosition(cellInfo.externalLineSegments(1:end, 1)), cellInfo.xPosition(cellInfo.externalLineSegments(1:end, 2))]',... 
       [cellInfo.yPosition(cellInfo.externalLineSegments(1:end, 1)), cellInfo.yPosition(cellInfo.externalLineSegments(1:end, 2))]',...
       'b');
  hold on
   
  plot([ cellInfo.xw(1:end-1) ; cellInfo.xw(2:end) ]',...
       [ cellInfo.yw(1:end-1); cellInfo.yw(2:end) ]',...
       'k');
  axis manual
  
  if isfield(cellInfo, "modelType") && cellInfo.modelType == "timeStepper"
    quiver(cellInfo.xPosition(:), cellInfo.yPosition(:),...
           .2*cellInfo.xVelocity(:), .2*cellInfo.yVelocity(:),...
           0, 'k');
  end
  
  hold off
end