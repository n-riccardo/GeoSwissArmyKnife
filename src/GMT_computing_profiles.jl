#using GMT
#using LibGEOS
#using CSV
#using DataFrames
#using DelimitedFiles

const global where_my_functions_are = @__DIR__

include(where_my_functions_are*"/RhumbLinesCalculations.jl")

"""
compute\\_GNSS\\_profiles:

*Last Update:*
06-02-2025

*Author:* 
Riccardo Nucci (riccardo.nucci4@unibo.it)

*Description:*
This function compute the parallel and orthogonal velocities of GNSS stations along a profile.
Outputs are ready for GMT plotting.

*Required arguments:*
- cross\\_profile\\_half\\_width::Real  # half width of the profile in km
- start\\_end\\_points::Vector{<:Real}=[0.0, 0.0, 0.0, 0.0], # format is lon1,lat1,lon2,lat2
- velo\\_data1::Matrix{<:Real}=Array{Float64}(undef, 0, 0) # matrix in format lon,lat,ve,vn,vu,se,sn,su

*Optional arguments:*
- outer\\_output\\_folder::String # where intermediate results will be saved
- string\\_id::String             # a string to identify the profile
- string\\_id\\_GNSS::String=""   # a string to identify the GNSS velocity dataset
- num\\_points_track::Int64=500   # number of points along the track on which crossprofiles will be erected
- StationNames::Vector{String}    # vector containing the names of the stations

*Notes:*
- The function work in spherical Earth approximation
- Correlation between GNSS velocity components is ignored

"""
function compute_GNSS_profiles(cross_profile_half_width::Real, start_end_points::Vector{<:Real}, velo_data1::Matrix{<:Real}; outer_output_folder::String="", string_id::String="", string_id_GNSS::String="", num_points_track::Int64=500, StationNames::Vector{String}=String[])

    if isempty(velo_data1)
        print("\n-> Nothing to compute: please specify velo_data1::Matrix{<:Real} in the format lon-lat-ve-vn-vu-se-sn-su \n")
        return
    end

    output_folder=""

    if(outer_output_folder!="")

        output_folder = outer_output_folder * "/" * "Profile_"*string_id * "/GNSS_"*string_id_GNSS*"/"
        print("\n-> Output files will be saved in:" * output_folder * "\n")
        
        if (isdir(output_folder))
            rm(output_folder, recursive=true)
            print("\n-> Previous folder " * output_folder * " has been deleted\n")
        end

        mkpath(output_folder)

    end
        
    total_distance, azimuth_, dataOnTrack, distances, lon_rect, lat_rect = compute_tracks(cross_profile_half_width,start_end_points,num_points_track)

    # DataFrames relative to the SwathProfile
    Info=DataFrame(Lon_Start=[start_end_points[1]], Lat_Start=[start_end_points[2]], Lon_End=[start_end_points[3]], Lat_End=[start_end_points[4]], Total_Distance=[total_distance], Azimuth=[azimuth_], Width=[cross_profile_half_width*2])
    LonLatCentralTrack=DataFrame(Distances=distances, Lon=dataOnTrack[:,1], Lat=dataOnTrack[:,2])
    LonLatRect=DataFrame(Lon=lon_rect, Lat=lat_rect)

    if(output_folder!="")
        CSV.write(output_folder*"Info.csv", Info, writeheader=true)
        CSV.write(output_folder*"LonLatCentralTrack.csv", LonLatCentralTrack, writeheader=true)
        CSV.write(output_folder*"LonLatSwathProfile.csv", LonLatRect, writeheader=true)
    end

    # Fasten the computation, first select stations using a crude lon-lat rectangle:
    lon_min=minimum(lon_rect);
    lon_max=maximum(lon_rect);
    lat_min=minimum(lat_rect);
    lat_max=maximum(lat_rect);

    velo_data=copy(velo_data1);
    
    cond=((velo_data[:,1] .>= lon_min) .& (velo_data[:,1] .<= lon_max)) .& ((velo_data[:,2] .>= lat_min) .& (velo_data[:,2] .<= lat_max));
    
    velo_data=velo_data[cond,:];
    if(!isempty(StationNames))
        StationNames=StationNames[cond]
    end

    # Keep stations within the polygon:
    lon_lat_rect=[lon_rect lat_rect]
    polycontour=[lon_lat_rect[k,:] for k in 1:size(lon_lat_rect,1)];
    SwathPolygon=LibGEOS.Polygon([polycontour]);
    withinPoly=Bool[]
    for i=1:length(velo_data[:,1])
        onePoint=LibGEOS.Point(velo_data[i,1],velo_data[i,2]);
        push!(withinPoly,LibGEOS.intersects(onePoint,SwathPolygon));
    end

    velo_data=velo_data[withinPoly,:];
    if(!isempty(StationNames))
        StationNames=StationNames[withinPoly]
    end

    # Now compute the distance of each station from the beginning of the SwathProfile and compute parallel and orthogonal velocities
    proj_points_lon=Float64[];
    proj_points_lat=Float64[];
    disances_velo=Float64[];
    
    for i=1:length(velo_data[:,1])
        proj_point_lon,proj_point_lat=rhxrh(velo_data[i,1], velo_data[i,2], azimuth_, start_end_points[1], start_end_points[2], azimuth_-90)
        push!(proj_points_lon,proj_point_lon)
        push!(proj_points_lat,proj_point_lat)
        push!(disances_velo,rhumb_distance([velo_data[i,1], velo_data[i,2]], [proj_point_lon, proj_point_lat]))
    end
    
    dim_data_r = length(velo_data[:,1])
    velocity_parallel_data = zeros(dim_data_r)
    sigma_parallel_data = zeros(dim_data_r)
    velocity_orthogonal_data = zeros(dim_data_r)
    sigma_orthogonal_data = zeros(dim_data_r)

    # to avoid confusion, explicitly define the variables
    ve_data=velo_data[:,3];
    vn_data=velo_data[:,4];
    vu_data=velo_data[:,5];
    se_data=velo_data[:,6];
    sn_data=velo_data[:,7];
    su_data=velo_data[:,8];

    for i in 1:dim_data_r
        velocity_parallel_data[i] = ve_data[i] * sind(azimuth_) + vn_data[i] * cosd(azimuth_)
        velocity_orthogonal_data[i] = ve_data[i] * cosd(azimuth_) - vn_data[i] * sind(azimuth_)
        sigma_parallel_data[i] = sqrt(sind(azimuth_)^2 * se_data[i]^2 + cosd(azimuth_)^2 * sn_data[i]^2)
        sigma_orthogonal_data[i] = sqrt(cosd(azimuth_)^2 * se_data[i]^2 + sind(azimuth_)^2 * sn_data[i]^2)
    end

    # DataFrames relative to the velocities
    if(isempty(StationNames))
        StationNames=string.(1:length(velo_data[:,1]))
    end
    Velocities=DataFrame(Station=StationNames, Lon=velo_data[:,1], Lat=velo_data[:,2], Distance=disances_velo,  Ve=ve_data, Vn=vn_data, Vu=vu_data, Se=se_data, Sn=sn_data, Su=su_data, Vel_Parallel=velocity_parallel_data, Sigma_Parallel=sigma_parallel_data, Vel_Orthogonal=velocity_orthogonal_data, Sigma_Orthogonal=sigma_orthogonal_data, Lon_Proj=proj_points_lon, Lat_Proj=proj_points_lat)

    # Sort the dataframe based on Distance:
    sort!(Velocities, :Distance)

    if(output_folder!="")
        CSV.write(output_folder*"Velocities.csv", Velocities, writeheader=true)
    end

    return Info,LonLatCentralTrack,LonLatRect,Velocities

end
### (: ###
### (: ###
### (: ###
### (: ###
"""
compute\\_profile\\_from\\_grid:

*Last Update:*
06-02-2025

*Author:* 
Riccardo Nucci (riccardo.nucci4@unibo.it)

*Description:*
This function uses GMT.jl functions to compute the scalar profile and crossprofile from a grid. 

*Required arguments:*
- cross\\_profile\\_half\\_width::Real 
- start\\_end_points::Vector{<:Real}=[0.0, 0.0, 0.0, 0.0], 
- grid\\_directory::String #a string of the path to the grid file in netcdf format

*Optional arguments:*
- outer\\_output\\_folder::String # where intermediate results will be saved
- string\\_id::String             # a string to identify the profile
- string\\_id\\_grid::String=""   # a string to identify the grid
- num\\_points\\_track::Int64=500
- sampling\\_dist::Real=1         # distance in km between points along the crossprofile

*Notes:*
- The function work in spherical Earth approximation
- Crossprofiles are erected using great circle lines
"""
function compute_profile_from_grid(cross_profile_half_width::Real, start_end_points::Vector{<:Real}=[0.0, 0.0, 0.0, 0.0],
    grid_directory::String=""; outer_output_folder::String="", string_id::String="", string_id_grid::String="",num_points_track::Int64=500, sampling_dist::Real=1)

    if isempty(grid_directory)
        print("\n-> Nothing to compute: please specify grid_directory::String file in netcdf format \n")
        return
    end

    output_folder=""

    if(outer_output_folder!="")

        output_folder = outer_output_folder * "/" * "Profile_"*string_id * "/Grid_"*string_id_grid*"/"
        print("\n-> Output files will be saved in:" * output_folder * "\n")
        
        if (isdir(output_folder))
            rm(output_folder, recursive=true)
            print("\n-> Previous folder " * output_folder * " has been deleted\n")
        end

        mkpath(output_folder)

    end

    total_distance, azimuth_, dataOnTrack, distances, lon_rect, lat_rect = compute_tracks(cross_profile_half_width,start_end_points,num_points_track)

    # DataFrames relative to the SwathProfile
    Info=DataFrame(Lon_Start=[start_end_points[1]], Lat_Start=[start_end_points[2]], Lon_End=[start_end_points[3]], Lat_End=[start_end_points[4]], Total_Distance=[total_distance], Azimuth=[azimuth_], Width=[cross_profile_half_width*2])
    LonLatCentralTrack=DataFrame(Distances=distances, Lon=dataOnTrack[:,1], Lat=dataOnTrack[:,2])
    LonLatRect=DataFrame(Lon=lon_rect, Lat=lat_rect)

    if(output_folder!="")
        CSV.write(output_folder*"Info.csv", Info, writeheader=true)
        CSV.write(output_folder*"LonLatCentralTrack.csv", LonLatCentralTrack, writeheader=true)
        CSV.write(output_folder*"LonLatSwathProfile.csv", LonLatRect, writeheader=true)
    end

    Lon_Lat_points=hcat(LonLatCentralTrack.Lon, LonLatCentralTrack.Lat);

    # Problem: if the point is outside the grid, it will not be included in the track
    GRD_track = grdtrack(Lon_Lat_points, grid=grid_directory)

    Lon_track = GRD_track[:, 1]
    Lat_track = GRD_track[:, 2]
    Values_track = GRD_track[:, 3]

    # Correction: look at the equal values of Lon-Lat
    Values_with_NaN = zeros(length(Lon_Lat_points[:, 1]))

    for i in 1:length(Lon_Lat_points[:, 1])

        my_index = findall(((Lon_track .== Lon_Lat_points[i, 1]) .& (Lat_track .== Lon_Lat_points[i, 2])))

        if (length(my_index) > 1)
            print("-> ERROR: multiple matches should not happen\n")
        elseif (length(my_index) == 1)
            Values_with_NaN[i] = Values_track[my_index[1]]
        else
            Values_with_NaN[i] = NaN
        end

    end

    # Save the results relative to the central track
    dfCentralTrack = DataFrame(
    Distances = LonLatCentralTrack.Distances,
    Longitude = Lon_Lat_points[:, 1],
    Latitude = Lon_Lat_points[:, 2],
    Values = Values_with_NaN
    )

    width_profile=cross_profile_half_width*2;
    cross_profile = grdtrack(Lon_Lat_points, grid=grid_directory, crossprofile="$(width_profile)k/$(sampling_dist)k")
    cross_profile_dataframe = DataFrame(cross_profile)

    dfCrossProfile=DataFrame(Longitude = Float64[], Latitude = Float64[], Values = Float64[])

    for i in 1:length(Lon_Lat_points[:,1])

        lon_temp = cross_profile_dataframe[i, 1]
        lat_temp = cross_profile_dataframe[i, 2]
        values_temp = cross_profile_dataframe[i, 5] #checked on GMT documentation

        for j in 1:length(lon_temp)
            push!(dfCrossProfile, hcat(lon_temp[j], lat_temp[j], values_temp[j]))
        end

        push!(dfCrossProfile, hcat(NaN,NaN,NaN))

    end

    if(output_folder!="")
        CSV.write(output_folder*"/ProfileCentralTrack.csv", dfCentralTrack, writeheader=true)
        CSV.write(output_folder*"/CrossProfileValues.csv", dfCrossProfile, writeheader=false)
    end

    return Info, LonLatCentralTrack, LonLatRect, dfCentralTrack, dfCrossProfile
    
end
### (: ###
### (: ###
### (: ###
### (: ###
"""
    compute\\_parorth\\_grid:

    *Author:*
    Riccardo Nucci (riccardo.nucci4@unibo.it)

    *Description:*
    It computes the parallel and orthogonal velocities for a given azimuth and return them

"""
function compute_parorth_grid(ve_model::Matrix{<:Real}, vn_model::Matrix{<:Real}, azimuth_::Real)


    # Velocity calculations
    velocity_parallel_model = ve_model .* sind.(azimuth_) .+ vn_model .* cosd.(azimuth_)
    velocity_orthogonal_model = ve_model .* cosd.(azimuth_) .- vn_model .* sind.(azimuth_)

    return velocity_parallel_model, velocity_orthogonal_model

end
### (: ###
### (: ###
### (: ###
### (: ###
"""
compute\\_seismicity\\_profiles:

*Last Update:*
06-02-2025

*Author:* 
Riccardo Nucci (riccardo.nucci4@unibo.it)

*Required arguments:*
- cross\\_profile\\_half\\_width::Real  # half width of the profile in km
- start\\_end\\_points::Vector{<:Real}=[0.0, 0.0, 0.0, 0.0], # format is lon1,lat1,lon2,lat2
- seismicity\\_data1::Matrix{<:Real}=Array{Float64}(undef, 0, 0) # matrix in format lon, lat, depth, Magnitude

*Optional arguments:*
- outer\\_output\\_folder::String # where intermediate results will be saved
- string\\_id::String             # a string to identify the profile
- string\\_id\\_seism::String=""   # a string to identify the seismicity dataset
- num\\_points_track::Int64=500   # number of points along the track on which crossprofiles will be erected

*Notes:*
- The function work in spherical Earth approximation
"""
function compute_seismicity_profiles(cross_profile_half_width::Real, start_end_points::Vector{<:Real}=[0.0, 0.0, 0.0, 0.0], seismicity_data1::Matrix{<:Real}=Array{Float64}(undef, 0, 0); outer_output_folder::String, string_id::String, string_id_seism::String="",num_points_track::Int64=500)

    if isempty(seismicity_data1)
        print("\n-> Nothing to compute: please specify seismicity::Matrix{<:Real} in the format lon-lat-depth-magnitude \n")
        return
    end

    output_folder=""

    if(outer_output_folder!="")

        output_folder = outer_output_folder * "/" * "Profile_"*string_id * "/Seism_"*string_id_seism*"/"
        print("\n-> Output files will be saved in:" * output_folder * "\n")
        
        if (isdir(output_folder))
            rm(output_folder, recursive=true)
            print("\n-> Previous folder " * output_folder * " has been deleted\n")
        end

        mkpath(output_folder)

    end

    total_distance, azimuth_, dataOnTrack, distances, lon_rect, lat_rect = compute_tracks(cross_profile_half_width,start_end_points,num_points_track)

    # DataFrames relative to the SwathProfile
    Info=DataFrame(Lon_Start=[start_end_points[1]], Lat_Start=[start_end_points[2]], Lon_End=[start_end_points[3]], Lat_End=[start_end_points[4]], Total_Distance=[total_distance], Azimuth=[azimuth_], Width=[cross_profile_half_width*2])
    LonLatCentralTrack=DataFrame(Distances=distances, Lon=dataOnTrack[:,1], Lat=dataOnTrack[:,2])
    LonLatRect=DataFrame(Lon=lon_rect, Lat=lat_rect)

    if(output_folder!="")
        CSV.write(output_folder*"Info.csv", Info, writeheader=true)
        CSV.write(output_folder*"LonLatCentralTrack.csv", LonLatCentralTrack, writeheader=true)
        CSV.write(output_folder*"LonLatSwathProfile.csv", LonLatRect, writeheader=true)
    end

    # Fasten the computation, first select events using a crude lon-lat rectangle:
    lon_min=minimum(lon_rect);
    lon_max=maximum(lon_rect);
    lat_min=minimum(lat_rect);
    lat_max=maximum(lat_rect);

    seismicity_data=copy(seismicity_data1);

    cond=((seismicity_data[:,1] .>= lon_min) .& (seismicity_data[:,1] .<= lon_max)) .& ((seismicity_data[:,2] .>= lat_min) .& (seismicity_data[:,2] .<= lat_max));
    
    seismicity_data=seismicity_data[cond,:];

    # Keep events within the polygon:
    lon_lat_rect=[lon_rect lat_rect]
    polycontour=[lon_lat_rect[k,:] for k in 1:size(lon_lat_rect,1)];
    SwathPolygon=LibGEOS.Polygon([polycontour]);
    withinPoly=Bool[]
    for i=1:length(seismicity_data[:,1])
        onePoint=LibGEOS.Point(seismicity_data[i,1],seismicity_data[i,2]);
        push!(withinPoly,LibGEOS.intersects(onePoint,SwathPolygon));
    end

    seismicity_data=seismicity_data[withinPoly,:];

    # Now compute the distance of each event from the beginning of the SwathProfile
    proj_points_lon=Float64[];
    proj_points_lat=Float64[];
    disances_seism=Float64[];
    
    for i=1:length(seismicity_data[:,1])
        proj_point_lon,proj_point_lat=rhxrh(seismicity_data[i,1], seismicity_data[i,2], azimuth_, start_end_points[1], start_end_points[2], azimuth_-90)
        push!(proj_points_lon,proj_point_lon)
        push!(proj_points_lat,proj_point_lat)
        push!(disances_seism,rhumb_distance([seismicity_data[i,1], seismicity_data[i,2]], [proj_point_lon, proj_point_lat]))
    end

    # DataFrames relative to the seismicity
    Seismicity=DataFrame(Lon=seismicity_data[:,1], Lat=seismicity_data[:,2], Distance=disances_seism,  Depth=seismicity_data[:,3], Magnitude=seismicity_data[:,4], Lon_Proj=proj_points_lon, Lat_Proj=proj_points_lat)

    # Sort the dataframe based on Distance:
    sort!(Seismicity, :Distance)

    if(output_folder!="")
        CSV.write(output_folder*"Seismicity.csv", Seismicity, writeheader=true)
    end

end
### (: ###
### (: ###
### (: ###
### (: ###
"""
percentile\\_from\\_cross\\_sections:

*Author:* 
Riccardo Nucci (riccardo.nucci4@unibo.it)

*Description:*
Automatically compute the desired percentile in orth directions from CrossProfileValues.

*Required arguments:*
- my\\_percentile::Real # The percentile to compute (one value only)

*Optional arguments (specify one of the two):*
- path\\_to\\_file::String       # Path of CrossProfileValues.csv
- CrossProfileValues::DataFrame  # DataFrame containing the crossprofile values

*Notes:*
- NaN values are ignored in the computation of the percentile
"""
function percentile_from_cross_sections(my_percentile::Real;path_to_file::String="",CrossProfileValues::DataFrame=DataFrame())

    if(isempty(path_to_file) && isempty(CrossProfileValues))
        print("\n-> Nothing to compute: please specify path_to_file::String or CrossProfileValues::DataFrame \n")
        return 0
    end

    if(isempty(path_to_file))
        path_to_file=where_my_functions_are*"/Temp/temp.csv"
        CSV.write(path_to_file,CrossProfileValues,writeheader=false)
    end

    # Read the file into an array
    data = readdlm(path_to_file, ',', Float64)  

    # Initialize storage for blocks
    blocks = []
    current_block = []

    # Separate data into blocks
    for row in eachrow(data)
        if all(isnan.(row))  # Check for NaN row as block separator
            if !isempty(current_block)
                push!(blocks, vcat(current_block...))  # Save current block
                current_block = []  # Reset for the next block
            end
        else
            push!(current_block, row')  # Add row to the current block
        end
    end

    # Add the last block if it's non-empty
    if !isempty(current_block)
        push!(blocks, vcat(current_block...))
    end

    # Compute the desired percentile for the third column in each block
    desired_percentile = my_percentile  
    percentile_results = Float64[]

    for block in blocks
        if !isempty(block)
            third_column = block[:, 3]  # Extract third column
			valid_values=third_column[.!(isnan.(third_column))]
			 if isempty(valid_values)
				push!(percentile_results, NaN)  # If no valid values, add NaN
			else
				push!(percentile_results, quantile(valid_values, desired_percentile / 100))
			end
        end
    end

    return percentile_results
end

function compute_tracks(HalfWidth,start_end_points,num_points_track)

    # General info on the profile:
    total_distance=rhumb_distance(start_end_points[1:2], start_end_points[3:4])
    azimuth_=Eval_azimuth(start_end_points[1:2], start_end_points[3:4])

    # Lon-lat coordinates of the "num_points_track" points on the profile:
    incremental_distance=total_distance/(num_points_track+1)
    dataOnTrackGMT=GMT.sample1d([start_end_points[1] start_end_points[2];start_end_points[3] start_end_points[4]], resample="R+l", inc=string(incremental_distance)*"k");
    dataOnTrack=Matrix(dataOnTrackGMT)

    # Distance of these points from the starting point:
    distances=zeros(length(dataOnTrack[:,1]))
    for i=1:length(dataOnTrack[:,1])
        distances[i]=rhumb_distance([start_end_points[1],start_end_points[2]], [dataOnTrack[i,1],dataOnTrack[i,2]])
    end

    # Lon-lat coordinates of the rectangle representing the Swath profile:
    lat_rect=Float64[];
    lon_rect=Float64[];
    push!(lon_rect,dataOnTrack[1,1]);
    push!(lat_rect,dataOnTrack[1,2]);
    for i=1:length(distances)

        lon_temp=dataOnTrack[i,1];
        lat_temp=dataOnTrack[i,2];

        # print([lon_temp,lat_temp])

        lon_rect_temp,lat_rect_temp=destination([lon_temp,lat_temp],azimuth_-90, HalfWidth);

        push!(lon_rect,lon_rect_temp);
        push!(lat_rect,lat_rect_temp);

    end
    for i=1:length(distances)

        lon_temp=dataOnTrack[end-i+1,1];
        lat_temp=dataOnTrack[end-i+1,2];

        lon_rect_temp,lat_rect_temp=destination([lon_temp,lat_temp],azimuth_+90, HalfWidth);

        push!(lon_rect,lon_rect_temp);
        push!(lat_rect,lat_rect_temp);

    end
    push!(lon_rect,dataOnTrack[1,1]);
    push!(lat_rect,dataOnTrack[1,2]);

    return total_distance,azimuth_,dataOnTrack,distances,lon_rect,lat_rect
end