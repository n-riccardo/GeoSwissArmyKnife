function LabelFocalMechanisms(Mrr,Mtt,Mff,Mrt,Mrf,Mtf)

fthrust=zeros(length(Mrr))
fstrikeslip=zeros(length(Mrr))
fnormal=zeros(length(Mrr))

for i=1:length(Mrr)

    mrr=Mrr[i]
    mtt=Mtt[i]
    mff=Mff[i]
    mrt=Mrt[i]
    mrf=Mrf[i]
    mtf=Mtf[i]
	
    # remove the isotropic part
    E=(mrr + mtt + mff)/3
    Mmatrix=[mrr-E mrt mrf;
             mrt mtt-E mtf;
             mrf mtf mff-E]
    # compute eigenvalues and eigenvectors
    F = eigen(Mmatrix)
    Indices=sortperm(F.values)
    EigenvaluesOrdered=F.values[Indices]
    EigenvectorsOrdered=F.vectors[:, Indices]

    MinEigenValue=EigenvaluesOrdered[1]
    InterEigenValue=EigenvaluesOrdered[2]
    MaxEigenValue=EigenvaluesOrdered[3]

    # Eigenvectors
    Pvector=EigenvectorsOrdered[:,1]
    Bvector=EigenvectorsOrdered[:,2]
    Tvector=EigenvectorsOrdered[:,3]

    #compute the angle of each eigenvector w.r.t. the horizontal plane (t-f)
    Pangle=atan(Pvector[1],sqrt(Pvector[2]^2+Pvector[3]^2))
    if(Pangle<0)
        Pangle=Pangle+pi
    end
    Interangle=atan(Bvector[1],sqrt(Bvector[2]^2+Bvector[3]^2))
    if(Interangle<0)
        Interangle=Interangle+pi
    end
    Tangle=atan(Tvector[1],sqrt(Tvector[2]^2+Tvector[3]^2))
    if(Tangle<0)
        Tangle=Tangle+pi
    end

    fthrust[i]=sin(Tangle)^2
    fstrikeslip[i]=sin(Interangle)^2
    fnormal[i]=sin(Pangle)^2

end

condThrust=(fthrust .> fstrikeslip) .& (fthrust .> fnormal);
condStrikeSlip=(fstrikeslip .> fthrust) .& (fstrikeslip .> fnormal);
condNormal=(fnormal .> fthrust) .& (fnormal .> fstrikeslip);

return condThrust, condStrikeSlip, condNormal

end