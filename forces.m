function force = forces(cellInfo)
  force = calcAllForces([cellInfo.xPosition, cellInfo.yPosition],cellInfo);
  sum(abs(force))
end