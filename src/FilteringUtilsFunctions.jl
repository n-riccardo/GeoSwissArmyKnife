# using CSV
# using DataFrames
# using Distances
# using GeographicLib
# using Ipopt
# using JuMP
# using StatsBase
# using Printf

"""
FilteringUtilsFunctions is a small set of useful tools for managing GPS velocity datasets.
"""
# Define the Earth's radius
# https://www.jpz.se/Html_filer/wgs_84.html
global const Earth_Radius = 6371.0087714 #km

"""
    WeightStations(myDataset::DataFrame, thresholdNnb::Real; dispInfo::Bool=true) -> DataFrame

# Author:
Riccardo Nucci (riccardo.nucci4@unibo.it)

# Description:
This function merges stations that are closer than a given threshold. The weighted mean of the velocities is computed, where the weights are the inverse of the velocity variance.  
The function returns the dataset with the merged stations, and the names of modified stations are indicated with an "M" at the end.

# Arguments:
- `myDataset::DataFrame`: The dataset containing station data. The format **must** be as follows:
    - `Long`: Longitude of the station  
    - `Lat`: Latitude of the station  
    - `E_Rate`: East velocity of the station  
    - `N_Rate`: North velocity of the station  
    - `U_Rate`: Up velocity of the station  
    - `σ_E`: East velocity uncertainty of the station  
    - `σ_N`: North velocity uncertainty of the station  
    - `σ_U`: Up velocity uncertainty of the station  
    - `ρ`: Correlation of East and North velocities (not used)  
    - `Site`: Name of the station  

- `thresholdNnb::Real`: The threshold distance (in km) used to determine whether two stations are considered neighbors.  

# Optional Arguments:
- `dispInfo::Bool`: If `true`, the function prints diagnostic information (default: `true`).  

# Returns:
- `DataFrame`: The dataset with merged stations.  

# Notes:
- Currently, the function does not consider the correlation between velocities.  
- Spherical Earth approx. is used

"""
function WeightStations(myDataset::DataFrame, thresholdNnb::Real; dispInfo::Bool=true)

    print("-> WeightStations function (29/11/24) \n")
    print("-> Your threshold is ", thresholdNnb, " km \n")
    
    dimInputData = size(myDataset, 1)
    dimInputDataWeight = size(myDataset, 1)

    # This is needed to properly compute the spatial mean of station positions:
    tempListWeightsLonLat = ones(dimInputData)

    i = 1
    while true
        j = 1
        while true
            if (haversine([myDataset.Long[i],myDataset.Lat[i]], [myDataset.Long[j],myDataset.Lat[j]],  Earth_Radius) <= thresholdNnb) && (i != j)

                w_vei = 1 / (myDataset.σ_E[i]^2)
                w_vej = 1 / (myDataset.σ_E[j]^2)
                vei = myDataset.E_Rate[i]
                vej = myDataset.E_Rate[j]

                ve = (vei * w_vei + vej * w_vej) / (w_vei + w_vej)
                se = 1 / sqrt(w_vei + w_vej)

                w_vni = 1 / (myDataset.σ_N[i]^2)
                w_vnj = 1 / (myDataset.σ_N[j]^2)
                vni = myDataset.N_Rate[i]
                vnj = myDataset.N_Rate[j]

                vn = (vni * w_vni + vnj * w_vnj) / (w_vni + w_vnj)
                sn = 1 / sqrt(w_vni + w_vnj)

                w_vui = 1 / (myDataset.σ_U[i])^2
                w_vuj = 1 / (myDataset.σ_U[j])^2
                vui = myDataset.U_Rate[i]
                vuj = myDataset.U_Rate[j]

                vu = (vui * w_vui + vuj * w_vuj) / (w_vui + w_vuj)
                su = 1 / sqrt(w_vui + w_vuj)

                w_lon_lat_i=tempListWeightsLonLat[i]
                w_lon_lat_j=tempListWeightsLonLat[j]

                mean_lon = (myDataset.Long[i]*w_lon_lat_i + myDataset.Long[j]*w_lon_lat_j) / (w_lon_lat_i+w_lon_lat_j)
                mean_lat = (myDataset.Lat[i]*w_lon_lat_i + myDataset.Lat[j]*w_lon_lat_j) / (w_lon_lat_i+w_lon_lat_j)

                if(dispInfo)
                    # Inform about the stations that have been weighted
                    print("###############################################\n")
                    print("Stations ", myDataset.Site[i], " and ", myDataset.Site[j], " have been weighted and merged\n")
                    print("New Station ", myDataset.Site[i][1:end-1]*"M", " has been created\n")
                    print("###############################################\n")
                end

                # Update the dataset; the name of modified stations is indicated with an "M" at the end of the name
                myDataset[i,:] = [mean_lon, mean_lat, ve, vn, vu, se, sn, su, 0.0, myDataset.Site[i][1:end-1]*"M"]
                deleteat!(myDataset, j)

                # Update the weights of lon and lat for the new station
                tempListWeightsLonLat[i]=w_lon_lat_i+w_lon_lat_j
                deleteat!(tempListWeightsLonLat, j)

                dimInputDataWeight = size(myDataset, 1)

                # Go one position back because the dimension of the dataset has been reduced of one row
                j = j - 1

                if(j<i) #it should not happen
                    i = i - 1
                end

            end

            j = j + 1

            if j > dimInputDataWeight
                break
            end

        end

        i = i + 1

        if i > dimInputDataWeight
            break
        end

    end

    # Print some diagnostic information
    print("Finished. Overview:\n")
    print("There were ", dimInputData, " stations in the dataset\n")
    print("Now there are ", dimInputDataWeight, " stations in the dataset\n")
    print("There are ", dimInputData - dimInputDataWeight, " stations weighted in the dataset")

    return myDataset

end

"""
    compute_distance_matrix(lon::Array{Float64,1}, lat::Array{Float64,1}) -> Matrix{Float64}

# Author:
Riccardo Nucci (riccardo.nucci4@unibo.it)

# Description:
Computes the distance between each pair of stations.

# Arguments:
- `lon::Array{Float64,1}`: Longitudes of the stations.  
- `lat::Array{Float64,1}`: Latitudes of the stations.  

# Returns:
- `distance_matrix::Matrix{Float64}`: A matrix containing the distances between stations.  

# Notes:
- Spherical Earth approx. is used

"""
function compute_distance_matrix(lon,lat)

	n_stations=length(lon)
	
	distance_matrix=zeros(n_stations,n_stations)

	for i in 1:n_stations, j in 1:n_stations

		distance_matrix[i,j]=haversine([lon[i],lat[i]], [lon[j],lat[j]],  Earth_Radius) #result in km

	end

	return distance_matrix

end

function compute_distance_matrix(lon1,lat1,lon2,lat2)

	n_stations1=length(lon1)
    n_stations2=length(lon2)
	
	distance_matrix=zeros(n_stations1,n_stations2)

	for i in 1:n_stations1, j in 1:n_stations2

		distance_matrix[i,j]=haversine([lon1[i],lat1[i]], [lon2[j],lat2[j]],  Earth_Radius) #result in km

	end

	return distance_matrix

end

"""
    km2deg(distance_km::Float64) -> Float64

# Author:
Riccardo Nucci (riccardo.nucci4@unibo.it)

# Description:
Convert a distance from kilometers to degrees of latitude in a Spherical Earth.

"""
function km2deg(distance_km::Float64)::Float64

    degrees_per_km = 360.0 / (2 * π * Earth_Radius)

    return distance_km * degrees_per_km

end

"""
    km2deg_longitude(distance_km::Float64) -> Float64

# Author:
Riccardo Nucci (riccardo.nucci4@unibo.it)

# Description:
Convert a distance from kilometers to degrees of longitude in a Spherical Earth at a certain latitude.

"""
function km2deg_longitude(distance_km::Float64, latitude::Float64)::Float64
    
    lat_rad = deg2rad(latitude)
    
    # Length of one degree of longitude at the given latitude
    km_per_deg_longitude = Earth_Radius * cos(lat_rad) * 2 * π / 360
    
    # Convert distance in kilometers to degrees of longitude
    deg_longitude = distance_km / km_per_deg_longitude
    
    return deg_longitude

end

"""
    find_variability_of_the_velocity_field(my_dataset::DataFrame, distance_of_investigation::Float64; perc_down::Int64=25, perc_up::Int64=75) 
    -> Tuple{Array{Float64,1}, Array{Float64,1}, Array{Float64,1}, Array{Int64,1}}

# Author:
Riccardo Nucci (riccardo.nucci4@unibo.it)

# Description:
Finds the variability of the velocity field within a characteristic length around each station.  
The variability is defined as the difference between the third and first quartiles of the velocities of stations closer than the **Lc**.

# Arguments:
- `my_dataset::DataFrame`: The dataset containing station data. The format **must** be as follows:
    - `Long`: Longitude of the station  
    - `Lat`: Latitude of the station  
    - `E_Rate`: East velocity of the station  
    - `N_Rate`: North velocity of the station  
    - `U_Rate`: Up velocity of the station  
    - `σ_E`: East velocity uncertainty of the station (not used)  
    - `σ_N`: North velocity uncertainty of the station (not used)  
    - `σ_U`: Up velocity uncertainty of the station (not used)  
    - `ρ`: Correlation of East and North velocities (not used)  
    - `Site`: Name of the station (not used)  

- `Lc::Float64`: The characteristic distance (in km) within which stations are considered for variability computation.  

# Optional Arguments:
- `perc_down::Int64`: The lower percentile (default: `25`).  
- `perc_up::Int64`: The upper percentile (default: `75`).  

# Returns:
- `Δve::Array{Float64,1}`: Bulk variability of the East velocities.  
- `Δvn::Array{Float64,1}`: Bulk variability of the North velocities.  
- `Δvu::Array{Float64,1}`: Bulk variability of the Up velocities.  
- `analyzed_stat::Array{Int64,1}`: Number of stations analyzed for each case.  
"""
function find_variability_of_the_velocity_field(my_dataset, Lc, perc_down=25, perc_up=75 ; honly=false)

    print("-> Find the variability of the velocity within a certain distance (29/11/24)\n")
    print("-> Your distance is ", Lc, " km\n")
	print("-> Your percentiles: ", perc_down, " ",perc_up,"\n")

    Δve=zeros(size(my_dataset, 1))
    Δvn=zeros(size(my_dataset, 1))
	if(!honly)
		Δvu=zeros(size(my_dataset, 1))
	else
		Δvu=NaN .* zeros(size(my_dataset, 1))
	end
    analyzed_stat=zeros(size(my_dataset, 1))

    for i in 1:size(my_dataset, 1)

        lon_i= my_dataset.Long[i]
        lat_i= my_dataset.Lat[i]

        indices_proximity=Int64[]

        # Find stations closer than Lc (including the station itself)
        for j in 1:size(my_dataset, 1)

            lon_j= my_dataset.Long[j]
            lat_j= my_dataset.Lat[j]

            d = haversine((lon_i, lat_i), (lon_j, lat_j), Earth_Radius)

            if(d<=Lc)        
                push!(indices_proximity, j)
            end

        end

        # Velocities of the stations closer than the Lc
        ve_proximity = my_dataset.E_Rate[indices_proximity]
        vn_proximity = my_dataset.N_Rate[indices_proximity]

        number_for_statistics=length(ve_proximity)

        first_quartile_ve=quantile(ve_proximity, perc_down/100)
        third_quartile_ve=quantile(ve_proximity, perc_up/100)

        first_quartile_vn=quantile(vn_proximity, perc_down/100)
        third_quartile_vn=quantile(vn_proximity, perc_up/100)


        Δve[i]=third_quartile_ve-first_quartile_ve
        Δvn[i]=third_quartile_vn-first_quartile_vn
		
		if(!honly)
			vu_proximity = my_dataset.U_Rate[indices_proximity]
			first_quartile_vu=quantile(vu_proximity, perc_down/100)
			third_quartile_vu=quantile(vu_proximity, perc_up/100)
			Δvu[i]=third_quartile_vu-first_quartile_vu
		end
		
        analyzed_stat[i]=number_for_statistics

    end

    print("Variability found\n")

    return Δve, Δvn, Δvu, analyzed_stat

end

"""
    smooth_the_field(my_dataset::DataFrame, Δv_dataset::DataFrame, L_c::Float64; coverage_factor::Int64=4, min_analyzed_stat::Int64=2) 
    -> Tuple{DataFrame, Array{Float64,1}, Array{Float64,1}, Array{Float64,1}}

# Author:
Riccardo Nucci (riccardo.nucci4@unibo.it)

# Description:
This function helps detect anomalous velocities in the dataset by solving an optimization problem.  
The objective is to find velocities of the stations that are most compatible with the bulk variability analyzed at a certain characteristic length within the uncertainties.

# Arguments:
- `my_dataset::DataFrame`: The dataset containing station data. The format **must** be as follows:
    - `Long`: Longitude of the station  
    - `Lat`: Latitude of the station  
    - `E_Rate`: East velocity of the station  
    - `N_Rate`: North velocity of the station  
    - `U_Rate`: Up velocity of the station  
    - `σ_E`: East velocity uncertainty of the station  
    - `σ_N`: North velocity uncertainty of the station  
    - `σ_U`: Up velocity uncertainty of the station  
    - `ρ`: Correlation of East and North velocities (not used)  
    - `Site`: Name of the station (not used)  

- `Δv_dataset::DataFrame`: The dataset containing the variability of the velocities. The format **must** be as follows:
    - `Δve`: Variability of the East velocities  
    - `Δvn`: Variability of the North velocities  
    - `Δvu`: Variability of the Up velocities  
    - `analyzed_stat`: The number of stations analyzed for each case  

- `L_c::Float64`: Distance (in km) used to calculate the variability of the velocities.  

# Optional Arguments:
- `coverage_factor::Int64`: The coverage factor (default: `4`).  
- `min_analyzed_stat::Int64`: The minimum number of stations for the analysis (default: `2`).  

# Returns:
- `my_smoothed_dataset::DataFrame`: The dataset containing the smoothed velocities.  
- `k_values_ve::Array{Float64,1}`: The k values for the East velocities.  
- `k_values_vn::Array{Float64,1}`: The k values for the North velocities.  
- `k_values_vu::Array{Float64,1}`: The k values for the Up velocities.  

"""
function smooth_the_field(my_dataset, Δv_dataset, Lc, coverage_factor=1.5, min_analysed_stat=4 ; honly=true)

    print("-> Smooth the velocity field (29/11/24)\n <-")
    print("Your Lc is ", Lc, " km\n")

    my_smoothed_dataset=copy(my_dataset)

    D=compute_distance_matrix(my_dataset.Long, my_dataset.Lat)

    μ_e=my_dataset.E_Rate;
    μ_n=my_dataset.N_Rate;
	if(!honly)
		μ_u=my_dataset.U_Rate;
		μ=[μ_e;μ_n;μ_u]
	else
		μ=[μ_e;μ_n]
	end

    Δve=Δv_dataset.Δve
    Δvn=Δv_dataset.Δvn
	if(!honly)
		Δvu=Δv_dataset.Δvu
		Δv=[Δve;Δvn;Δvu]
	else
		Δv=[Δve;Δvn]
	end

    σ_e=my_dataset.σ_E;
    σ_n=my_dataset.σ_N;
	if(!honly)
		σ_u=my_dataset.σ_U;
		σ=[σ_e;σ_n;σ_u]
	else
		σ=[σ_e;σ_n]
	end

    analyzed_stat=Δv_dataset.analyzed_stat

    #number_of_stations
    ns=Int(length(μ_e))

    model = Model(optimizer_with_attributes(Ipopt.Optimizer))

    # Define variables
	if(!honly)
		@variable(model, x[1:(3*ns)]) # ve, vn, vu
	else
		@variable(model, x[1:(2*ns)]) # ve, vn, vu
	end

    # Define objective function: minimize the sum of the squared differences between the smoothed velocities and the original velocities
    @objective(model, Min, sum(((x-μ).^2.0)./σ.^2.0))

    k_values_ve=Float64[]
    k_values_vn=Float64[]
    k_values_vu=Float64[]

    # Add nonlinear constraints
    for i=1:ns
        for j=1:ns

            if(j<i)

                if(D[i,j]<Lc)

                    if((analyzed_stat[i]>=min_analysed_stat) & (analyzed_stat[j]>=min_analysed_stat))

                        k_e_i=abs(coverage_factor*Δv[i])/Lc                
                        k_n_i=abs(coverage_factor*Δv[i+ns])/Lc
                        
                        k_e_j=abs(coverage_factor*Δv[j])/Lc
                        k_n_j=abs(coverage_factor*Δv[j+ns])/Lc
                        
                        k_e=max(k_e_i,k_e_j)
                        k_n=max(k_n_i,k_n_j)
                        
                        k_values_ve=push!(k_values_ve,k_e)
                        k_values_vn=push!(k_values_vn,k_n)
					
                    
                        @constraint(model, 
                        ((x[i]-x[j])^2)/(D[i,j]^2) <= (k_e^2))

                        @constraint(model, 
                        ((x[ns+i]-x[ns+j])^2)/(D[i,j]^2) <= (k_n^2))

						if(!honly)
							k_u_i=abs(coverage_factor*Δv[i+2*ns])/Lc
							k_u_j=abs(coverage_factor*Δv[j+2*ns])/Lc
							k_u=max(k_u_i,k_u_j)
							k_values_vu=push!(k_values_vu,k_u)
							@constraint(model, 
							((x[2*ns+i]-x[2*ns+j])^2)/(D[i,j]^2) <= (k_u^2))
						end

                    end
                end
            end
        end
    end

    # Solve the optimization problem
    optimize!(model)

    smoothed_velocities=value.(x)
    ve_smoothed=smoothed_velocities[1:ns]
    vn_smoothed=smoothed_velocities[(ns+1):(2*ns)]

    my_smoothed_dataset.E_Rate=ve_smoothed
    my_smoothed_dataset.N_Rate=vn_smoothed
	
	if(!honly)
		vu_smoothed=smoothed_velocities[(2*ns+1):end]
		my_smoothed_dataset.U_Rate=vu_smoothed
	end

    return my_smoothed_dataset, k_values_ve, k_values_vn, k_values_vu

end

"""
    find_suspect_outliers(my_dataset::DataFrame, my_smoothed_datset::DataFrame, use_uncertainties::Bool=false) -> Array{Float64,1}

# Author:
Riccardo Nucci (riccardo.nucci4@unibo.it)

# Description:
This function compute the compatibility of the original velocities with respect to the smoothed velocities. It return a `log_likelihood` value for each station, representing the estimate of compatibility.

# Arguments:
- `my_dataset::DataFrame`: The dataset containing station data. The format **must** be as follows:
    - `Long`: Longitude of the station  
    - `Lat`: Latitude of the station  
    - `E_Rate`: East velocity of the station  
    - `N_Rate`: North velocity of the station  
    - `U_Rate`: Up velocity of the station  
    - `σ_E`: East velocity uncertainty of the station  
    - `σ_N`: North velocity uncertainty of the station  
    - `σ_U`: Up velocity uncertainty of the station  
    - `ρ`: Correlation of East and North velocities (not used)  
    - `Site`: Name of the station (not used)  
- `my_smoothed_datset::DataFrame`: The dataset containing the smoothed velocity field. The format **must** be as follows:
    - `E_Rate`: East velocity of the station  
    - `N_Rate`: North velocity of the station  
    - `U_Rate`: Up velocity of the station  

# Optional Arguments:
- `use_uncertainties::Bool`: If `true`, the function uses the uncertainties of the velocities to compute the likelihood (default: `false`).

"""
function find_suspect_outliers(my_dataset,my_smoothed_datset,use_uncertainties=false)
    
    ve=my_dataset.E_Rate
    vn=my_dataset.N_Rate
    #vu=my_dataset.U_Rate

    sigma_e=my_dataset.σ_E
    sigma_n=my_dataset.σ_N
    #sigma_u=my_dataset.σ_U

    ve_smoothed=my_smoothed_datset.E_Rate
    vn_smoothed=my_smoothed_datset.N_Rate
    #vu_smoothed=my_smoothed_datset.U_Rate

    log_likelihoods=zeros(Float64,length(ve))

    for i in 1:length(my_dataset[:,1])

        ve_temp=ve[i]
        vn_temp=vn[i]
        #vu_temp=vu[i]

        ve_smoothed_temp=ve_smoothed[i]
        vn_smoothed_temp=vn_smoothed[i]
        #vu_smoothed_temp=vu_smoothed[i]

        sigma_e_temp=sigma_e[i]
        sigma_n_temp=sigma_n[i]
        #sigma_u_temp=sigma_u[i]

        if(use_uncertainties)
            likelihood=exp(-0.5*(((ve_smoothed_temp-ve_temp)/sigma_e_temp)^2+((vn_smoothed_temp-vn_temp)/sigma_n_temp)^2))
        else
            likelihood=exp(-0.5*(((ve_smoothed_temp-ve_temp))^2+((vn_smoothed_temp-vn_temp))^2))
        end

        log_likelihoods[i]=log10(likelihood)

    end

    return log_likelihoods

end

"""
    RemoveByCircularArea(my_dataset::DataFrame, coords_and_radius::Matrix{Float64}) -> FinalDataset::DataFrame

*Author*:
Riccardo Nucci (riccardo.nucci4@unibo.it)

*Description*:
Stations within a certain radius from a given point are removed from the dataset. The function returns the dataset without the stations in the circular area.

*Required arguments:*
- my\\_dataset::DataFrame: the dataset containing the stations. The format MUST be the following:
    - Long: Longitude of the station
    - Lat: Latitude of the station
    - E_Rate: East velocity of the station
    - N_Rate: North velocity of the station
    - U_rate: Up velocity of the station
    - σ_E: East velocity uncertainty of the station
    - σ_N: North velocity uncertainty of the station
    - σ_U: Up velocity uncertainty of the station
    - ρ: Correlation of East and North velocities (not used)
    - Site: Name of the station (not used)
- coords\\_and\\_radius::Matrix{Float64}: a matrix containing in each row: longitude, latitude, and radius of the circular area
"""
function RemoveByCircularArea(myDataset, coords_and_radius)

    dim_of_input_data_v = size(myDataset, 1)
    indices_to_be_cleaned = Int[]

    for j in 1:size(coords_and_radius, 1)
        for i in 1:dim_of_input_data_v

            temp_distance=haversine([myDataset.Long[i],myDataset.Lat[i]], [coords_and_radius[j,1],coords_and_radius[j,2]],  Earth_Radius) 
            
            if temp_distance <= coords_and_radius[j,3]
                push!(indices_to_be_cleaned, i)
            end

        end
    end

    # Delete selected rows
    FinalDataset = myDataset[Not(indices_to_be_cleaned), :]

    return FinalDataset,indices_to_be_cleaned

end

function generate_circles(coords_and_radius; num_points=360)
    circles = []

    for j in 1:size(coords_and_radius, 1)

        lon_center = coords_and_radius[j, 1]
        lat_center = coords_and_radius[j, 2]
        radius_km = coords_and_radius[j, 3]

        lattrk = Float64[]
        lontrk = Float64[]

        for azimuth in range(0, 360, length=num_points)
            lon, lat, _ = GeographicLib.forward(lon_center, lat_center, azimuth, radius_km * 1000)
            push!(lontrk, lon)
            push!(lattrk, lat)
        end

        push!(circles, (lontrk, lattrk))  # Store each circle as a tuple (lon_list, lat_list)
    end

    return circles
end

function RemoveByName(my_dataset, namelist)

    dim_of_input_data_v = size(my_dataset, 1)
    indices_to_be_cleaned = Int[]
    list_of_found_stat=String[]

    for j in 1:length(namelist)
        for i in 1:dim_of_input_data_v

            if my_dataset.Site[i] == namelist[j]
                push!(indices_to_be_cleaned, i)
                push!(list_of_found_stat, namelist[j])
            end

        end
    end

    print("Number of found station(s): "*string(length(indices_to_be_cleaned)) * " out of "*string(length(namelist))*"\n")
    print("I do not find the following station(s): "*string(setdiff(namelist,list_of_found_stat))*"\n")
    my_dataset = my_dataset[Not(indices_to_be_cleaned), :]

    return my_dataset

end

function RemoveByPosition(my_dataset, coords, tolerance_lon, tolerance_lat)

    dim_of_input_data_v = size(my_dataset, 1)
    indices_to_be_cleaned = Int[]

    for j in 1:size(coords, 1)
        for i in 1:dim_of_input_data_v

            if abs(my_dataset.Long[i] - coords[j,1]) <= tolerance_lon && abs(my_dataset.Lat[i] - coords[j,2]) <= tolerance_lat
                push!(indices_to_be_cleaned, i)
            end
        end
    end
    
    my_dataset = my_dataset[Not(indices_to_be_cleaned), :]

    return my_dataset

end

"""
# Arguments:
- MasterDataset, Dataframe, Fomat: lon,lat,ve,vn,vu,se,sn,su,name
- SlaveDataset, Dataframe, Fomat: lon,lat,ve,vn,vu,se,sn,su,name
- Dtreshold: distance threshold in km
"""
function RotateVeloField(MasterDataset, SlaveDataset, Dtreshold; components=3)

    # Find the common sites
    distanceMatrix=compute_distance_matrix(MasterDataset[:,1],MasterDataset[:,2],SlaveDataset[:,1],SlaveDataset[:,2])
    cond=(distanceMatrix .< Dtreshold);

    # find rows and columns of the common sites
    couple_of_indices=Tuple.(findall(!iszero,cond))

    MasterIndices= first.(couple_of_indices)
    SlaveIndices= last.(couple_of_indices)

    # Construct the design matrix
    δV=Float64[]
    R=Matrix{Float64}(undef,0,3)
    O=Float64[]

    for i in eachindex(view(cond, :, 1))
        for j in eachindex(view(cond, 1, :))
            if cond[i,j]

                lon=MasterDataset[i,1]
                lat=MasterDataset[i,2]
                δvE=MasterDataset[i,3]-SlaveDataset[j,3]
                δvN=MasterDataset[i,4]-SlaveDataset[j,4]

                RE=Earth_Radius* [sind(lon) -cosd(lon) 0];
                RN=Earth_Radius* [-cosd(lon)*sind(lat) -sind(lon)*sind(lat) cosd(lat)];

                if(components==3)
                    δvU=MasterDataset[i,5]-SlaveDataset[j,5]
                    RU=[0 0 0]
                    Oij=[0,0,1]

                    R=vcat(R,RE,RN,RU)

                    push!(δV,δvE)
                    push!(δV,δvN)
                    push!(δV,δvU)

                    append!(O,Oij)
                else
                    R=vcat(R,RE,RN)
                    push!(δV,δvE)
                    push!(δV,δvN)
                end

            end
        end
    end

    if components==3
        G=hcat(R,O)
    else
        G=R
    end
    
    # Get the parameters through least square solution
    w = G \ δV
    
    # Now rotote Slave into Master
    vEnew=Float64[]
    vNnew=Float64[]
    vUnew=Float64[]
    vErot=Float64[]
    vNrot=Float64[]
    vUrot=Float64[]
    
    for i in eachindex(view(SlaveDataset, :, 1))
        
        lon=SlaveDataset[i,1]
        lat=SlaveDataset[i,2]

        vEold=SlaveDataset[i,3]
        vNold=SlaveDataset[i,4]
    
        RE=Earth_Radius* [sind(lon) -cosd(lon) 0];
        RN=Earth_Radius* [-cosd(lon)*sind(lat) -sind(lon)*sind(lat) cosd(lat)];
    
        if(components==3)

            vUold=SlaveDataset[i,5]
            RU=[0 0 0]
            O=[0,0,1]
    
            R=vcat(RE,RN,RU)
            Gtemp=hcat(R,O)

            vnewtemp=Gtemp*w + [vEold,vNold,vUold]
            push!(vEnew,vnewtemp[1])
            push!(vNnew,vnewtemp[2])
            push!(vUnew,vnewtemp[3])

            vrottemp=Gtemp*w;
            push!(vErot,vrottemp[1])
            push!(vNrot,vrottemp[2])
            push!(vUrot,vrottemp[3])

        else

            Gtemp=vcat(RE,RN)

            vnewtemp=Gtemp*w + [vEold,vNold]
            push!(vEnew,vnewtemp[1])
            push!(vNnew,vnewtemp[2])

            vrottemp=Gtemp*w;
            push!(vErot,vrottemp[1])
            push!(vNrot,vrottemp[2])

        end
    
    end

    if components==3
        SlaveDataset.vERotated=vEnew
        SlaveDataset.vnRotated=vNnew
        SlaveDataset.vuRotated=vUnew

        SlaveDataset.RotVE=vErot
        SlaveDataset.RotVN=vNrot
        SlaveDataset.RotVU=vUrot
    else
        SlaveDataset.vERotated=vEnew
        SlaveDataset.vnRotated=vNnew

        SlaveDataset.RotVE=vErot
        SlaveDataset.RotVN=vNrot
    end

    return MasterDataset, SlaveDataset, MasterIndices, SlaveIndices

end

function call_velrot(MasterDataset, SlaveDataset, Dtreshold; param_opt="TR",vert_weight=1)

    # parm_opt: TRS is translation, rotation and scaling
    
    SlaveFile=where_my_functions_are*"/../temp/Slave.vel"
    MasterFile=where_my_functions_are*"/../temp/Master.vel"

    write_eq_dist_velrot(where_my_functions_are*"/../temp/link.file", Dtreshold)

    # Write Master in velrot format
    open(MasterFile, "w") do io
        for i in 1:length(MasterDataset[:,1])
            My_Vector = [MasterDataset[i,1], MasterDataset[i,2], 
                         MasterDataset[i,3], MasterDataset[i,4], 
                         0.0, 0.0, MasterDataset[i,6], MasterDataset[i,7], 
                         0.0, MasterDataset[i,5], 0.0, MasterDataset[i,8], 
                         MasterDataset[i,9]]
            println(io, format_row_velrot(My_Vector))
        end
    end

    # Write Slave in velrot format
    open(SlaveFile, "w") do io
        for i in 1:length(SlaveDataset[:,1])
            My_Vector = [SlaveDataset[i,1], SlaveDataset[i,2], 
                         SlaveDataset[i,3], SlaveDataset[i,4], 
                         0.0, 0.0, SlaveDataset[i,6], SlaveDataset[i,7], 
                         0.0, SlaveDataset[i,5], 0.0, SlaveDataset[i,8], 
                         SlaveDataset[i,9]]
            println(io, format_row_velrot(My_Vector))
        end
    end

    # Call velrot that rotates sys1 (slave) into sys2 (master)
    cd(where_my_functions_are*"/../temp/") do 
        run(`velrot Slave.vel eura Master.vel eura output_velrot.txt eura link.file $vert_weight $param_opt`)
    end

    outVelrotName=where_my_functions_are*"/../temp/output_velrot.txt";

    # Reading section
    # Read the output file
    indx_fsites, indx_SlaveInMaster, indx_MasterInSlave = get_indices(outVelrotName)

    # Residuals:
    ResDataset = CSV.read(outVelrotName, 
    header=["*", "#", "Name_1", "Name_Ref", "dN", "dE", "dU", "sN", "sE", "sU", "sTN", "sTE", "sTU"], 
    skipto=indx_fsites + 2, limit=indx_SlaveInMaster - indx_fsites - 3, delim=' ', ignorerepeated=true, DataFrame)

    # Common stations
    indxCommonSlave = zeros(Bool, length(SlaveDataset[:,1]))
    indxCommonMaster = zeros(Bool, length(MasterDataset[:,1]))

    CommonSlaveRep=DataFrame()
    CommonMasterRep=DataFrame()

    for i in 1:length(ResDataset.Name_1)

        name_temp_Slave = ResDataset.Name_1[i]
        name_temp_Master = ResDataset.Name_Ref[i]

        for j in 1:length(SlaveDataset[:,1])
            if SlaveDataset[j,9] == name_temp_Slave
                indxCommonSlave[j]=true
                push!(CommonSlaveRep, SlaveDataset[j,:])
            end
        end
        for k in 1:length(MasterDataset[:,1])
            if MasterDataset[k,9] == name_temp_Master
                indxCommonMaster[k]=true
                push!(CommonMasterRep, MasterDataset[k,:])
            end
        end
    end

    #CommonSlave=SlaveDataset[indxCommonSlave.==true,:]
    #CommonMaster=MasterDataset[indxCommonMaster.==true,:]

    # Rotated Datasets
    SlaveInMaster = CSV.read(outVelrotName,
     header=["Long", "Lat", "E_Rate", "N_Rate", "E_adj", "N_adj", "σ_E", "σ_N", "ρ", "U_Rate", "U_adj", "σ_U", "Site"],
      skipto=indx_SlaveInMaster + 3, limit=indx_MasterInSlave - indx_SlaveInMaster - 4, delim=' ', ignorerepeated=true, DataFrame)

    for row in 1:nrow(SlaveInMaster)
        # Remove the last character if it's '*' or '+'
        SlaveInMaster[row, :Site] = remove_last_char(SlaveInMaster[row, :Site])
    end

    temp_file = where_my_functions_are*"/../temp/temp_velrot.txt"

    copy_and_modify_file(outVelrotName, temp_file, indx_MasterInSlave)

    MasterInSlave = CSV.read(temp_file, 
    header=["Long", "Lat", "E_Rate", "N_Rate", "E_adj", "N_adj", "σ_E", "σ_N", "ρ", "U_Rate", "U_adj", "σ_U", "Site"],
    skipto=4, delim=' ', ignorerepeated=true, DataFrame)

    rm(temp_file)

    for row in 1:nrow(MasterInSlave)
        # Remove the last character if it's '*' or '+'
        MasterInSlave[row, :Site] = remove_last_char(MasterInSlave[row, :Site])
    end

    print(indx_MasterInSlave)
    return indxCommonMaster, indxCommonSlave, SlaveInMaster, MasterInSlave
end


function format_row_velrot(My_Vector)
    return @sprintf("%8.3f %8.3f   %6.2f   %6.2f   %6.2f   %6.2f   %5.2f   %5.2f  %6.3f    %6.2f  %6.2f  %6.2f  %s",
        round(My_Vector[1], digits=3), round(My_Vector[2], digits=3), round(My_Vector[3], digits=2), round(My_Vector[4], digits=2), 
        round(My_Vector[5], digits=2), round(My_Vector[6], digits=2), round(My_Vector[7], digits=2), round(My_Vector[8], digits=2), 
        round(My_Vector[9], digits=3), round(My_Vector[10], digits=2), round(My_Vector[11], digits=2), round(My_Vector[12], digits=2), 
        My_Vector[13])
end

function write_eq_dist_velrot(filename::String, distance::Real)
    distance=distance*1000.0;
    open(filename, "w") do io
        println(io, "# Set the nominal separation of sites to be this distance (m)")
        println(io, " eq_dist $distance")
    end
end

function remove_last_char(str)
    if endswith(str, '*') || endswith(str, '+')
        return chop(str)
    else
        return str
    end
end

function copy_and_modify_file(input_file::String, output_file::String, line_start::Int64)
    # Equal station names are indicated with "-". We don't need them.
    input_stream = open(input_file, "r")
    output_stream = open(output_file, "w")

    line_number = 1

    for line in eachline(input_stream)

        if line_number >= line_start
            modified_line = line[2:end]
            println(output_stream, modified_line)
        end
        line_number += 1
    end
    close(input_stream)
    close(output_stream)
end

function get_indices(output_velrot)
    
    file = open(output_velrot, "r")
    
    indx_fsites = -1
    indx_SlaveInMaster = -1
    indx_MasterInSlave = -1
    my_index = 1
    
    for line in eachline(file)
    
        if (contains(line, "*  Differences at the fundamental sites"))
            indx_fsites = my_index
        end
        if (contains(line, "* SYSTEM 1 Velocities transformed to SYSTEM 2 "))
            indx_SlaveInMaster = my_index
        end
        if (contains(line, "* SYSTEM 2 Velocities except those in SYSTEM 1 "))
            indx_MasterInSlave = my_index
        end
    
        my_index = my_index + 1
    end
    
    close(file)

    return indx_fsites, indx_SlaveInMaster, indx_MasterInSlave

end

function combine_velo_field(MasterDataset, SlaveDataset, Dtreshold)

    # Identical sites below a threshold:
    distanceMatrix=compute_distance_matrix(MasterDataset[:,1],MasterDataset[:,2],SlaveDataset[:,1],SlaveDataset[:,2])
    cond=(distanceMatrix .< Dtreshold);

    # find rows and columns of the common sites
    couple_of_indices=Tuple.(findall(!iszero,cond))

    #MasterIndices= first.(couple_of_indices)
    SlaveIndices= last.(couple_of_indices)

    # remove stations in the SlaveDataset that are in the MasterDataset
    SlaveDataset = SlaveDataset[Not(SlaveIndices), :]

    rename!(MasterDataset, [:lon, :lat, :ve, :vn, :vu, :se, :sn, :su, :name]);
    rename!(SlaveDataset, [:lon, :lat, :ve, :vn, :vu, :se, :sn, :su, :name]);
    
    # merge the two datasets
    FinalDataset = vcat(MasterDataset, SlaveDataset)

    # order by the last column
    FinalDataset = sort(FinalDataset,  order(:name))

    return FinalDataset

end
