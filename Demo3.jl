using StaticArrays
using LinearAlgebra
using Random
using GLMakie
using Purses
include("lib.jl")
include("cached.jl")

sel_n = 1
sel_j = 1
zeros = zeros(Float64,length(360))
contact_x = getindex.(contact_pts[1],1)
contact_y = getindex.(contact_pts[1],2)
contact_z = getindex.(contact_pts[1],3)
dev_x = real.(conn_dev[1])
dev_y = imag.(conn_dev[1])
sel_S = refS[1][sel_n]
sel_C = 0.8*sel_S
sel_circ = circ[1][sel_n]
sel_holom2 = HolomorphicTransform(conn_dev[1],sel_n)
sel_holom3 = -InvStereoProj.(sel_holom2)
holom2_x = Observable.(real.(sel_holom2))
holom2_y = Observable.(imag.(sel_holom2))
holom3_x = Observable.(getindex.(sel_holom3,1))
holom3_y = Observable.(getindex.(sel_holom3,2))
holom3_z = Observable.(getindex.(sel_holom3,3))
circ_x = Observable.(getindex.(sel_circ,1))
circ_y = Observable.(getindex.(sel_circ,2))
circ_z = Observable.(getindex.(sel_circ,3))

fig = Figure(size = (1600, 850))
Label(fig[1, 1],"Contact curve",fontsize = 20)
Label(fig[1, 2],"Development",fontsize = 20)
Label(fig[1, 3],"Holomorph2",fontsize = 20)
Label(fig[1, 4],"Holomorph3",fontsize = 20)

ax_contact=Axis3(fig[2, 1])
ax_dev=Axis3(fig[2, 2])
ax_holom2=Axis3(fig[2, 3])
ax_holom3=Axis3(fig[2, 4])

lines!(ax_contact,contact_x,contact_y,contact_z,linewidth = 2)
lines!(ax_dev,dev_x,dev_y,zeros,linewidth = 2)
lines!(ax_holom2,holom2_x,holom2_y,zeros,linewidth = 2)
lines!(ax_holom3,holom3_x,holom3_y,holom3_z,linewidth = 2)
lines!(ax_holom3,circ_x,circ_y,circ_z,linewidth = 1)


# ============================================================
# Slider
# ============================================================
n_slider = Slider(fig[4, 2:3],range = 1:360,startvalue = 1)
Label(fig[3, 2:3],lift(n -> "n = $n",n_slider.value),fontsize = 18)
# ============================================================
# Selected point
# ============================================================
scatter!(ax_dev,lift(n -> Point3f(dev_x[n],dev_y[n],zeros[n]),n_slider.value),color = :red,markersize = 15)
# ============================================================
# Slider callback
# ============================================================
on(n_slider.value) do n
  sel_n = n
  sel_holom2[] = HolomorphicTransform(conn_dev[1],sel_n)
  sel_holom3[] = -InvStereoProj.(sel_holom2)
  holom2_x[] = real.(sel_holom2)
  holom2_y[] = imag.(sel_holom2)
  holom3_x[] = getindex.(sel_holom3,1)
  holom3_y[] = getindex.(sel_holom3,2)
  holom3_z[] = getindex.(sel_holom3,3)
  circ_x[] = getindex.(sel_circ,1)
  circ_y[] = getindex.(sel_circ,2)
  circ_z[] = getindex.(sel_circ,3)

end

display(fig)






