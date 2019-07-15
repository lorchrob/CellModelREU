%{
Wrapper functino for plotCell (i.e., its only job is to call plotCell).
Set up to be called by the ode solvers.
%}
function status = plotCellWrapper(~, y, flag, externalNodeCount, simulationType, modelType)
  persistent n
  if isempty(n)
      n = 0;
  end
    
  % only plot every 15 steps  
  if mod(n, 15) == 0
    if strcmp(flag, 'init') || strcmp(flag, 'done')
      % do nothing
    else
      if exist('simulationType', 'var')
        cellInfo = initializeNetwork(externalNodeCount, simulationType);
      else
        cellInfo = initializeNetwork(externalNodeCount); 
      end
    
      if exist('modelType', 'var')
        cellInfo.modelType = modelType;
      end
    
      cellInfo.xPosition = y(1:2:end);
      cellInfo.yPosition = y(2:2:end);
      plotCell(cellInfo)
      drawnow;
    end
  end
  
  n = n + 1;
  status = 0;
end