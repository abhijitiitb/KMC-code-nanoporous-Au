%removeoverlap.m: removes solvent with overlap

%create linked list for solvent
ncellsLinkedList=ceil(BoxSize(1:3)/CellSizeLinkedList); %ecah cell is 4 Angstrom in size

%create LinkedListMetal
LinkedListMetal.N=zeros(ncellsLinkedList);
LinkedListMetal.particles=zeros([ncellsLinkedList 100]); %assume 100 particles in each cell
for iatom=1:NMetalAtoms
    icell=ceil(metalcoord(iatom,:)/CellSizeLinkedList);
    LinkedListMetal.N(icell(1),icell(2),icell(3))=LinkedListMetal.N(icell(1),icell(2),icell(3))+1;
    n=LinkedListMetal.N(icell(1),icell(2),icell(3));
    LinkedListMetal.particles(icell(1),icell(2),icell(3),n)=iatom;
end

%check for overlapping solvent
for iatom=1:NSolvent
    nooverlap=1;
    coords=solventcoord(iatom,:);
    icell=ceil(coords/CellSizeLinkedList);
    icell=max(icell,1);
    if any(icell==1) || any(icell==ncellsLinkedList) %no point is checking this particle
        %continue %as it is at the edge
    else
        for i=-1:1
            for j=-1:1
                for k=-1:1
                    jcell=icell+[i j k];
                    nooverlap=nooverlap & nooverlapfound(coords,metalcoord,LinkedListMetal,jcell,tolerance);
                end
            end
        end
        mask(iatom)=nooverlap;
    end
end
