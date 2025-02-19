"""
Script to create Hovmoller diagram of the the TCBL flow from Scythe output.

Tyler Barbero, February 2025
"""

using Statistics, CairoMakie, CSV, DataFrames

# define time and radius for Hovmoller
times = collect(0:3600:86400*4)
nt = length(times)
nr = 150 # length of radius array

# allocate memory for array
radius = Array{Float64}(undef,nr)
wbts_g = Array{Float64}(undef,nr,nt);
ubts_g = similar(wbts_g)
vbts_g = similar(wbts_g)
wbts_mr = similar(wbts_g)
ubts_mr = similar(wbts_g)
vbts_mr = similar(wbts_g)

# compile data into an array from each hourly file
for (i,timet) in enumerate(times)
    expname = "/bell-scratch/tbarbero/Scythe_test/R1D_test/yes_f/physical_out_"
    loaddir = expname*map(string,timet)*".0.csv"

    # open file and store data
    df = CSV.read(loaddir,DataFrame)
    wb = df.wb
    ub = df.ub
    vb = df.vb
    
    wbts_g[:,i] = wb
    ubts_g[:,i] = ub
    vbts_g[:,i] = vb
    
    expname = "/bell-scratch/tbarbero/Scythe_test/Williams2013_tcbl/modified_Rankine_0.3/physical_out_"
    loaddir = expname*map(string,timet)*".0.csv"
    
    # open file and store data
    df = CSV.read(loaddir,DataFrame)
    wb = df.wb
    ub = df.ub
    vb = df.vb
    
    wbts_mr[:,i] = wb
    ubts_mr[:,i] = ub
    vbts_mr[:,i] = vb
    
    if i==1
         global radius = df.r/1000.
    end
end

# Functions for plotting
#============================================================#
function create_axis(f_loc)
    xtl=["0","","200","","400","","600"]
    ax = Axis(f_loc,
        yticks=0:24:96,
        xticks=((0:100:600),xtl),
#         xticklabels=xtl,
        aspect=1.5,
    )
end

function plot_labels(ax,title="",xlabel="",ylabel="")
    xlims!(ax,0,600)
    ylims!(ax,0,96)
    ax.titlefont=:regular
    ax.title = title
    ax.ylabel = ylabel
    ax.xlabel = xlabel
end

function plot_ub(ax,x,y,z,clevels)
    return dat = contourf!(ax,x,y,z,
                    levels=clevels,
                    colormap=Reverse("deep"),
                    extendlow=:auto,
                    extendhigh=:auto
                    )
end

function plot_vb(ax,x,y,z,clevels)
    return dat = contourf!(ax,x,y,z,
                    levels=clevels,
                    colormap=:viridis,
                    extendhigh=:auto
                    )
end

function plot_wb(ax,x,y,z,clevels)
    dat = contourf!(ax,x,y,z,
                    levels=clevels,
                    colormap=:bam,
                    extendhigh=:auto
                    )
        
    # zero contours
    contour!(ax,x,y,z,
                    levels=[0],
                    color=:black,
                    linewidth=2)

    return dat
end

function plot_colorbar(f,dat,clevels)
    return Colorbar(f,
            dat,
            vertical=false,
            ticks=clevels)
end
#============================================================#

# main plotting
f = Figure(size=(1500,850),
            fontsize=25,
            figure_padding=35)

#============================#
ax = create_axis(f[1,1])
clevels_ub = -5:1:0
dat = plot_ub(ax,radius,times/3600.,ubts_g,clevels_ub)
plot_labels(ax,L"$u_b$","","Gaussian \n t(hr)")
plot_colorbar(f[3,1],dat,clevels_ub)

#============================#
ax = create_axis(f[1,2])
clevels_vb = 0:2:15
dat = plot_vb(ax,radius,times/3600.,vbts_g,clevels_vb)
# plot_labels(ax,title=L"$v_b$")
plot_labels(ax,L"$v_b$")
plot_colorbar(f[3,2],dat,clevels_vb)

#============================#
ax = create_axis(f[1,3])
clevels_wb = -10:2:10
dat = plot_wb(ax,radius,times/3600.,wbts_g*100,clevels_wb)
plot_labels(ax,L"100 * $w_b$")
plot_colorbar(f[3,3],dat,clevels_wb)

#============================#
ax = create_axis(f[2,1])
plot_ub(ax,radius,times/3600.,ubts_mr,clevels_ub)
plot_labels(ax,"","r(km)","modified Rankine \n t(hr)")

#============================#
ax = create_axis(f[2,2])
plot_vb(ax,radius,times/3600.,vbts_mr,clevels_vb)
plot_labels(ax,"","r(km)","")

#============================#
ax = create_axis(f[2,3])
plot_wb(ax,radius,times/3600.,wbts_mr*100,clevels_wb)
plot_labels(ax,"","r(km)","")

# save figure
savepath = "/home/tbarbero/Hovmoller_idealized_vortices_TCBL.png"
save(savepath,f)
