%{
Function to plot a cell with a red membrane and a blue interior
%}
function plotCell(cellInfo)  

  extLineSegs = cellInfo.externalLineSegments;
  lineSegs = cellInfo.lineSegments;
  IorE = logical(ismember(lineSegs, extLineSegs, 'rows') + ismember(flip(lineSegs,2), extLineSegs, 'rows'));
  
  plot( [cellInfo.xPosition(cellInfo.lineSegments(IorE, 1)), cellInfo.xPosition(cellInfo.lineSegments(IorE, 2))]',... 
       [cellInfo.yPosition(cellInfo.lineSegments(IorE, 1)), cellInfo.yPosition(cellInfo.lineSegments(IorE, 2))]',...
       'r');
  axis([-20 20 -20 20])
  hold on
  plot([cellInfo.xPosition(cellInfo.lineSegments(~IorE, 1)), cellInfo.xPosition(cellInfo.lineSegments(~IorE, 2))]',... 
       [cellInfo.yPosition(cellInfo.lineSegments(~IorE, 1)), cellInfo.yPosition(cellInfo.lineSegments(~IorE, 2))]',...
       'b');
  axis manual
  
  if isfield(cellInfo, "modelType") && cellInfo.modelType == "timeStepper"
     
    quiver(cellInfo.xPosition(:), cellInfo.yPosition(:),...
           cellInfo.xVelocity(:), cellInfo.yVelocity(:),...
           0, 'k');
  end
  
  hold off
end