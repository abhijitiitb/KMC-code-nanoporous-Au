%run.m: exclude solvent overlapping with the nanoporous metal structure
%it is assumed that the solvent box is larger than the metal, however, the
%two overlap

l=10; %Angstrom, solvent box is l A additional width in +x and -x compared to metal, similarly for y and z directions
CellSizeLinkedList=4; %typically l > CellSizeLinkedList ... see removeoverlap.m
tolerance=3; %Angstrom for removing overlap

%metal details
filemetal='2red.coordinates';
readmetal; %read metal file
BoxSizeMetal(4:6)=BoxSizeMetal(4:6)-l;  %New addition to script

for i=1:NMetalAtoms           %New addition to script
  metalcoord(i,1:3) = metalcoord(i,1:3) - BoxSizeMetal(4:6);     %New addition to script
end      %New addition to script
%metalcoord=metalcoord-BoxSizeMetal(4:6)+l; %shift metal appropriately so that solvent start at origin

BoxSizeMetal(4:6)=BoxSizeMetal(4:6)+l;  %New addition to script
BoxSizeMetal=[BoxSizeMetal(1:3)-BoxSizeMetal(4:6)+l l l l];

%solvent details
iscreatesolvent=1; %whether solvent need to be created
if iscreatesolvent %create solvent using cubic unit cell
    
    asolvent=3.1; %lattice parameter of solvent in Angstrom (based on water density)
    
    nunitcells=[BoxSizeMetal(1)-BoxSizeMetal(4)+2*l BoxSizeMetal(2)-BoxSizeMetal(5)+2*l BoxSizeMetal(3)-BoxSizeMetal(6)+2*l];
    nunitcells=ceil(nunitcells/asolvent);
    createsolvent; %note: solvent starts at origin
    
    
else %read file
    filesolvent='solvent.xyz'; %solvent file
    BoxSize=[10 10 10 0 0 0]; %specify size of the solvent box
    readsolvent; %reads the solvent file
    
end
mask=ones(NSolvent,1); %for solvent particles 1- retain solvent, 0- remove solvent

%write initial file
fileout='metal-solvent-overlap.xyz';
writeXYZ;

%remove overlap
removeoverlap;

%write final file
fileout='metal-solvent-without-overlap.xyz';
writeXYZ;

exit
