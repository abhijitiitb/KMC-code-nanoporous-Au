%readmetal.m: read coordinates in metal file

A=dlmread(filemetal);
NMetalAtoms=size(A,1);
metalcoord=A(:,3:5);
BoxSizeMetal=[max(metalcoord(:,1)) max(metalcoord(:,2)) max(metalcoord(:,3)) ...
    min(metalcoord(:,1)) min(metalcoord(:,2)) min(metalcoord(:,3)) ];
MetalSpecies=A(:,2);


