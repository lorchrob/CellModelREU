%{
Function to plot a cell, using a breadth-first search algorithm
%}
function plotCell(cellInfo)  
  for i = 1 : numel(cellInfo.lineSegments) / 2
     plot([cellInfo.xPosition(cellInfo.lineSegments(i, 1)), cellInfo.xPosition(cellInfo.lineSegments(i, 2))],... 
          [cellInfo.yPosition(cellInfo.lineSegments(i, 1)), cellInfo.yPosition(cellInfo.lineSegments(i, 2))],...
          'b');
      
     hold on
  end
  
  hold off
end