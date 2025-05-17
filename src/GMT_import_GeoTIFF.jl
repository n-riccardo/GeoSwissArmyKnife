# From Geophysical Model Generator:

function lonlatdepth_grid(Lon::Any, Lat::Any, Depth::Any)

  nLon    = length(Lon)
  nLat    = length(Lat)
  nDepth  = length(Depth)

  if nLon==nLat==nDepth==1
      error("Cannot use this routine for a 3D point (no need to create a grid in that case, 1 point only")
  end
  if maximum([length(size(Lon)), length(size(Lat)), length(size(Depth))])>1
      error("You can only give 1D vectors or numbers as input")
  end

  Lon3D   =   zeros(nLon,nLat,nDepth);
  Lat3D   =   zeros(nLon,nLat,nDepth);
  Depth3D =   zeros(nLon,nLat,nDepth);

  for i=1:nLon
      for j=1:nLat
          for k=1:nDepth
              Lon3D[i,j,k]    =   ustrip.(Lon[i]);
              Lat3D[i,j,k]    =   ustrip.(Lat[j]);
              Depth3D[i,j,k]  =   ustrip.(Depth[k]);
          end
      end
  end

  # Add dimensions back
  Lon3D   =   Lon3D*unit(  Lon[1])
  Lat3D   =   Lat3D*unit(  Lat[1])
  Depth3D = Depth3D*unit(Depth[1])

  return Lon3D, Lat3D, Depth3D
end

function import_GeoTIFF(fname::String; fieldname=:layer1, negative=false, iskm=true, NorthernHemisphere=true, constantDepth=false, removeNaN_z=false, removeNaN_field=false)
  G = gmtread(fname);

  # Transfer to GeoData
  nx,ny = length(G.x)-1, length(G.y)-1
  Lon,Lat,Depth   =   lonlatdepth_grid(G.x[1:nx],G.y[1:ny],0);
  if  hasfield(typeof(G),:z) 
    Depth[:,:,1]    =   G.z';
    if negative
      Depth[:,:,1]  =   -G.z';
    end
    if iskm
      Depth    *=   1e-3*km;
    end
  end

  # Create GeoData structure
  data = zero(Lon)
  if hasfield(typeof(G),:z)
    data = Depth
  
  elseif hasfield(typeof(G),:image)
    if length(size(G.image)) == 3
      data = permutedims(G.image,[2, 1, 3]);
    elseif length(size(G.image)) == 2
      data[:,:,1] = G.image'
    end

  end
  
  if removeNaN_z
    remove_NaN_surface!(Depth, Lon, Lat)
  end
  if removeNaN_field
    remove_NaN_surface!(data, Lon, Lat)
  end
  data_field  = NamedTuple{(fieldname,)}((data,));

  if constantDepth
    Depth = zero(Lon)
  end

  if contains(G.proj4,"utm")
    zone = parse(Int64,split.(split(G.proj4,"zone=")[2]," ")[1]); # retrieve UTM zone
    data_GMT    = UTMData(Lon, Lat, Depth, zone, NorthernHemisphere, data_field)
  
  elseif contains(G.proj4,"longlat") 
    data_GMT    = GeoData(Lon, Lat, Depth, data_field)

  else
    error("I'm sorry, I don't know how to handle this projection yet: $(G.proj4)\n
           We recommend that you transfer your GeoTIFF to longlat by using QGIS \n
           Open the GeoTIFF there and Export -> Save As , while selecting \"EPSG:4326 - WGS 84\" projection.")
  end

  return data_GMT
end