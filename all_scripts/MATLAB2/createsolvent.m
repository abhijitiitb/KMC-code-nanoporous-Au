%createsolvent.m: creates a cube of solvent molecules
solventcoord=zeros(prod(nunitcells),3);
NSolvent=0;
for icell=1:nunitcells(1)
    x=(icell-1)*asolvent;
    for jcell=1:nunitcells(2)
        y=(jcell-1)*asolvent;
        for kcell=1:nunitcells(3)
            z=(kcell-1)*asolvent;
            NSolvent=NSolvent+1;
            solventcoord(NSolvent,:)=[x y z];
        end
    end
end
BoxSize=[nunitcells*asolvent 0 0 0];