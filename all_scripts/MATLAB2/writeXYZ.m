%WriteXYZ.m: write out the xyz coordinates

fid=fopen(fileout,'w');
fprintf(fid,'%d\n',NMetalAtoms+sum(mask)); %sum(mask) gives the number of solvent
fprintf(fid,'%s\n','Overlapping particles');
for iatom=1:NMetalAtoms
    if MetalSpecies(iatom,1)==1
        species='Ag';
    else
        species='Au';
    end
    fprintf(fid,'%s %f %f %f\n',species,metalcoord(iatom,:));
end

for iatom=1:NSolvent
    species='W';
    if mask(iatom)
        fprintf(fid,'%s %f %f %f\n',species,solventcoord(iatom,:));
    end
end

fclose(fid);