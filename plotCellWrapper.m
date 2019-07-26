%{
Wrapper functino for plotCell (i.e., its only job is to call plotCell).
Set up to be called by the ode solvers.
%}
function status = plotCellWrapper(t, y, flag, externalNodeCount, simulationType, modelType)
  % 'n' used to keep track of the iteration so that we don't have to plot
  % every time
  
  % weird variable name intended to make sure the global variable
  % isn't affected by anything else
  persistent n12
  if isempty(n12)
      n12 = -1;
  end
  
  % aren't defined until the end, used for calculating velocities
  persistent yOld
  persistent tOld
  
  % only plot every 15 steps  
  if mod(n12, 5) == 0
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
    
      if exist('yOld', 'var') 
        cellInfo.xVelocity = (y(1:2:end,1) - yOld(1:2:end,1))/(t(1) - tOld(1));
        cellInfo.yVelocity = (y(2:2:end,1) - yOld(2:2:end,1))/(t(1) - tOld(1));
      end
      
      cellInfo.xPosition = y(1:2:end,1);
      cellInfo.yPosition = y(2:2:end,1);
      plotCell(cellInfo)
      drawnow;
    end
  end
  
  n12 = n12 + 1;
  status = 0;
  
  yOld = y;

  tOld = t;
end