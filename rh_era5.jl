#-----------------------------------------------------------------------------------------------
# rh_era5.jl
#
# script producing a single panel map of relative humidity (RH) between +/-40 degrees latitude.
# the RH that is plotted is the difference between El Nino and La Nina years. 
#
# the figures is saved as a file named: fig2name (a png file)
# 
# ERA5 data is used after having been interpolated to a 1x1 degree grid
# ask levi if you want the data file.  
#
# ensoFuncs.jl contains many of the functions used in this script and needs to be run before
# this script will work.
#
# defined outside of this script: 
# numfields, high, low
#
# the high and low arrays are created in CMIProniComposites.jl by calling the check_thresh_high
# and check_thresh_low functions.  if the lat/lon indices defined in CMIProniComposites.jl are
# not correct then the figures produced with this script will be wrong.   
# 
# levi silvers                                                                feb 2025
#-----------------------------------------------------------------------------------------------

using CairoMakie
using GeoMakie
using NCDatasets
using Statistics


# select the hemisphere to focus on.   these indices are based on the peak TC seasons
# for each hemisphere. 
# NH
#ninoyears = [18 54 90 150 174 234 306 402]
#ninayears = [102 114 210 246 318 366 378 390]
# SH
ninoyears = [23 35 59 96 155 239 311 347]
ninayears = [107 119 215 251 263 335 371 383]

fig2name = tag*"_rh_nino_comp_SH_test.png"

function create_indices(years)
  ensoInd = Matrix{Int64}(undef, 8, 6)
  for i in 1:8
    si = years[i]
    ensoInd[i,:] = si:si+5
  end
  return ensoInd
end 

ninoInd = create_indices(ninoyears)
ninaInd = create_indices(ninayears)
#

path="/Users/C823281551/data/ERA5/"

filein  = path*"era5_rh_1990th2023_360x80.nc"
tag = "ERA5"
data   = NCDataset(filein)

lat = data["lat"]
lon = data["lon"]
lev = data["pressure_level"]

rh_var = data["r"]

numfields = 48 
numall    = 408

lat1      = 1
lat2      = 81

dims = size(rh_var)

rh_1            = Array{Union{Missing, Float64}, 4}(undef, dims[1], dims[2], 1, numfields)
rh_2            = Array{Union{Missing, Float64}, 4}(undef, dims[1], dims[2], 1, numfields)
rh_low          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
rh_high         = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numfields)
rh_tot          = Array{Union{Missing, Float64}, 3}(undef, dims[1], dims[2], numall)

endi = 48
# low is an array that contains the timesteps representing the negative phase of ENSO
# high is an array that contains the timesteps representing the positive phase of ENSO
low = ninaInd
high = ninoInd

# select the case determining if rh will be plotted at one level or at the average of two
pcase = 0# pcase = 0 corresponds to plotting rh on 1 level
# calculate total rh field, for all times
level=2 # level 2 should correspond to the 700 hPa pressure level. 
if pcase == 0
  for i in 1:408
    rh_tot[:,lat1:lat2,i]   = rh_var[:,lat1:lat2,level,i] 
  end
  for i in 1:endi
    # low values, la nina
    rh_2[:,lat1:lat2,1,i]   = rh_var[:,lat1:lat2,level,low[i]]
    rh_low[:,lat1:lat2,i]   = rh_2[:,lat1:lat2,1,i] 
    # high values, el nino
    rh_2[:,lat1:lat2,1,i]   = rh_var[:,lat1:lat2,level,high[i]]
    rh_high[:,lat1:lat2,i]  = rh_2[:,lat1:lat2,1,i] 
  end
else
  for i in 1:408
    rh_tot[:,lat1:lat2,i]   = (rh_var[:,lat1:lat2,2,i] + rh_var[:,lat1:lat2,1,i])./2
  end
  for i in 1:endi
    # low values, la nina
    rh_1[:,lat1:lat2,1,i]   = rh_var[:,lat1:lat2,1,low[i]]
    rh_2[:,lat1:lat2,1,i]   = rh_var[:,lat1:lat2,2,low[i]]
    rh_low[:,lat1:lat2,i]   = (rh_2[:,lat1:lat2,1,i] + rh_1[:,lat1:lat2,1,i])./2
    # high values, el nino
    rh_1[:,lat1:lat2,1,i]   = rh_var[:,lat1:lat2,1,high[i]]
    rh_2[:,lat1:lat2,1,i]   = rh_var[:,lat1:lat2,2,high[i]]
    rh_high[:,lat1:lat2,i]  = (rh_2[:,lat1:lat2,1,i] + rh_1[:,lat1:lat2,1,i])./2
  end
end 


rh_tot_tmn  = mean(rh_tot, dims=3)
rh_high_tmn = mean(rh_high, dims=3)
rh_low_tmn  = mean(rh_low, dims=3)

data_2_plot_tot  = rh_tot_tmn
data_2_plot_anom = rh_high_tmn - rh_low_tmn

function fig_anom_plot(inpv,d1,d2,tit)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:white,
        #size=(500,200),
        #size=(800,300),
        size=(600,300),
        )
    #ax = Axis(f2[1,1]; #--> default plot is rectangular equidistant 
    ax = GeoAxis(f2[1,1];
        xticks = -180:30:180, 
        #xticks = 0:30:360, 
        yticks = -90:30:90,
        ylabel="latitude",
        xlabel="longitude",
        limits=(-180,180,-40,40),
        title=tit,
        )
        bb = contourf!(ax, d1, d2, inpv, 
             levels = range(-14, 14, length = 21), # rh
             colormap = :BrBg,
             extendlow = :auto, extendhigh = :auto
        )
        lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.75)
        Colorbar(f2[1,2], bb)
    return f2
end
function fig_tot_plot(inpv,d1,d2,tit)
    f2 = Figure(;
        figure_padding=(5,5,10,10),
        backgroundcolor=:white,
        size=(600,300),
        )
    ax = GeoAxis(f2[1,1];
        xticks = -180:30:180, 
        #xticks = 0:30:360, 
        yticks = -90:30:90,
        ylabel="latitude",
        xlabel="longitude",
        limits=(-180,180,-40,40),
        title=tit,
        )
        bb = contourf!(ax, d1, d2, inpv, 
             levels = range(0, 50, length = 20), # rh
             colormap = :bam,
             extendlow = :auto, extendhigh = :auto
        )
        lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.75)
        Colorbar(f2[1,2], bb)
    return f2
end

tit="ERA5 RH Composite (%)"
fig = fig_anom_plot(data_2_plot_anom[:,:,1],lon,lat,tit)

save(fig2name, fig)
