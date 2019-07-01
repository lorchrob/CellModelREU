%{
Function to plot a cell with a red membrane and a blue interior
%}
function plotCell(cellInfo)  
  for i = 1 : numel(cellInfo.lineSegments) / 2
     if cellInfo.lineSegments(i, 1) <= cellInfo.externalNodeCount && cellInfo.lineSegments(i, 2) <= cellInfo.externalNodeCount
       plot([cellInfo.xPosition(cellInfo.lineSegments(i, 1)), cellInfo.xPosition(cellInfo.lineSegments(i, 2))],... 
            [cellInfo.yPosition(cellInfo.lineSegments(i, 1)), cellInfo.yPosition(cellInfo.lineSegments(i, 2))],...
            'r');
     else
       plot([cellInfo.xPosition(cellInfo.lineSegments(i, 1)), cellInfo.xPosition(cellInfo.lineSegments(i, 2))],... 
            [cellInfo.yPosition(cellInfo.lineSegments(i, 1)), cellInfo.yPosition(cellInfo.lineSegments(i, 2))],...
            'b');    
     end
     hold on
  end
  hold off
end