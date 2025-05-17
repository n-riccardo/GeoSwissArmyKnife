function ComputeHorizontalStrainRate(Exx::Vector{<:Real}, Eyy::Vector{<:Real}, Exy::Vector{<:Real})

    MaxEigenValues=[]
    MinEigenValues=[]
    Azimuths=[]

    for i=1:length(Exx)

        if(isnan(Exx[i]))
            push!(MaxEigenValues,NaN)
            push!(MinEigenValues,NaN)
            push!(Azimuths,NaN)
        else
            EMatrixTemp=[Exx[i] Exy[i]; Exy[i] Eyy[i]]
            # compute eigenvalues and eigenvectors
            F = eigen(EMatrixTemp)

            Indices=sortperm(F.values)
            EigenvaluesOrdered=F.values[Indices]
            EigenvectorsOrdered=F.vectors[:, Indices]

            MinEigenValue=EigenvaluesOrdered[1]
            MaxEigenValue=EigenvaluesOrdered[2]

            # Compute the azimuth direction of the most compressional eigenvalue
            MostComprVector=EigenvectorsOrdered[:,1]

            azimuth=atan(MostComprVector[1],MostComprVector[2]);

            if(azimuth<0)
                azimuth=azimuth+2*pi
            end

            azimuthDeg=azimuth*180/pi

            push!(MaxEigenValues,MaxEigenValue)
            push!(MinEigenValues,MinEigenValue)
            push!(Azimuths,azimuthDeg)
        end

    end

    return MaxEigenValues, MinEigenValues, Azimuths

end

function ReadVisrOutput()
	
	
	
end


