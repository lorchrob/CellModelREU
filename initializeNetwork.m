%{ 
TO DO:
  * make structs to store information about the network
      - number of internal and external nodes
      - position   
      - physical properties such as reference lengths, pressure,
        viscoelasticity, etc.
  * generate mesh
%}

function [cellInfo, systemInfo, programInfo] = initializeNetwork(nodeCount)  
  cellInfo.nodeCount = nodeCount
end