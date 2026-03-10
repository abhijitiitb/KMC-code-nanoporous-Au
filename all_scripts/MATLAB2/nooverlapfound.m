function f=nooverlapfound(coords,metalcoord,LinkedListMetal,jcell,tolerance)

f=1;
for iatom=1:LinkedListMetal.N(jcell(1),jcell(2),jcell(3))
    matom=LinkedListMetal.particles(jcell(1),jcell(2),jcell(3),iatom);
    coordm=metalcoord(matom,:);
    v=coords-coordm;
    v=sqrt(dot(v,v));
    if v<tolerance
        f=0;
    end
end