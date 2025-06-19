#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot the points where RI first occurs on a global map between +/-40 degrees lat
#
# this script was used to create panel Figure 2, panel a in 
# Silvers and Klotzbach et al., 2025, GRL
#
# see various projection options here: 
# https://proj.org/en/9.4/operations/projections/index.html
#
# run with: 
# /Users/silvers/.juliaup/bin/julia tcRIscatterGRLFig.jl
#
# levi silvers                                                                      june 2025
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# /Users/silvers/data/hurricane_data/AtlanticTCACE18512023.csv 

using DataFrames
using CSV
using CairoMakie
using LaTeXStrings
using GeoMakie

# incoming file: 
filein = "/Users/C823281551/data/hurricane_data/FirstRILocationbyStorm.csv"
file2in = "/Users/C823281551/data/hurricane_data/RI_Episodes.csv"

# read a csv file into a dataframe:
#df1 = DataFrame(CSV.File(filein))

dfRI  = CSV.read(filein, DataFrame)
dfRIe = CSV.read(file2in, DataFrame)

# how do I look at the metadata for df?
describe(dfRI)
# one can also describe  particular columns: 
#describe(dfRI, cols=3:10)

#img = rotr90(GeoMakie.earth())

function fig_pan(xScRI,yScRI,intens,tit)
    figpan = Figure(;
             figure_padding=(5,5,10,10),
             backgroundcolor=:white,
             size=(600,300),
             )
    ax = GeoAxis(figpan[1,1];
        xticks = -180:30:180, 
        yticks = -90:30:90,
        ylabel="latitude",
        xlabel="longitude",
        limits=(-180,180,-40,40),
        title=tit,
        )
        #bb = scatter!(ax, yScRI,xScRI, 
        #      markersize = 4, 
        #      color = :black     # was heat 
        #)
        bb = scatter!(ax, yScRI,xScRI,
              color=intens,
              markersize = 4, # originally was 7
              colormap = :matter     # was heat 
        )
        #
        #bb = contourf!(ax, d1, d2, inpv, 
        #     levels = range(-20, 20, length = 20), # rh
        #     #colormap = :Blues_8,
        #     #colormap = :broc,
        #     colormap = :bam,
        #     #colormap = :batlow,
        #     #colormap = :vik,
        #     extendlow = :auto, extendhigh = :auto
        #)
        lines!(ax, GeoMakie.coastlines(), color = :black, linewidth=0.75)
        Colorbar(figpan[2,1], vertical = false, flipaxis = false, 
                 bb)
    return figpan
end
fig = Figure(
    size = (1920, 880),
)
#ga = GeoAxis(
#    fig[1, 1]; # any cell of the figure's layout
#    dest = "+proj=wintri", # the CRS in which you want to plot
#)
#ga = GeoAxis(
#image!(ga, -180..180, -90..90, img; interpolate = false) 
#lines!(ga, GeoMakie.coastlines()) # plot coastlines from Natural Earth as a reference

# looking at these storms (763-1110) corresponds roughly to the events between 2014-2014.
p1=1

## location of first RI
#p2=1110
#xRI = dfRI[p1:p2, :LATITUDE]
#yRI = dfRI[p1:p2, :LONGITUDE]

## all RI events
p2=4490
xRI = dfRIe[p1:p2, :USA_LAT]
yRI = dfRIe[p1:p2, :USA_LON]
zRI = dfRIe[p1:p2, :"24-Hr Intensification Rate"]
#yearRI = dfRI[:, :YEAR]

##scatter(yRI,xRI)
#scatter!(ga, yRI,xRI, markersize = 7, color = :black)
#
#ylims!(ga, -40, 40)
#
## select by basin: 
#dfNA=dfRI[in.(dfRI.BASIN, Ref(["NA"])), :]
#xRI = dfNA[:, :LATITUDE]
#yRI = dfNA[:, :LONGITUDE]
#yearRI = dfNA[:, :YEAR]
#scatter!(ga, yRI,xRI, markersize = 7, colormap = :heat)

tit = " "
fig = fig_pan(xRI,yRI,zRI,tit) 

#save("tcRIscatterGRLpanFirstRI.png", fig)
save("tcRIscatterGRLpanRI.png", fig)

