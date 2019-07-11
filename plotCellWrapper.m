%{
Wrapper functino for plotCell (i.e., its only job is to call plotCell).
Set up to be called by the ode solvers.
%}
function status = plotCellWrapper(~, y, flag, externalNodeCount)
  if strcmp(flag, 'init') || strcmp(flag, 'done')
  else
    cellInfo = initializeNetwork(externalNodeCount); 
    cellInfo.xPosition = y(1:2:end)';
    cellInfo.yPosition = y(2:2:end)';
    plotCell(cellInfo)
    drawnow;
  end
  status = 0;
end