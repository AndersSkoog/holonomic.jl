using StaticArrays
using LinearAlgebra
using Random
using GLMakie
include("Constructors.jl")

fig = Figure(size = (800,800))
#display contact sphere curve
ax_contact = Axis3(fig[4,1])
#display development curve
ax_dev = Axis3(fig[4,2])
#display transformed development curve with respect to n
ax_trans_dev = Axis3(fig[4,3])
#display the transformed development lifted to the sphere
ax_holomorph = Axis3(fig[4,4])

n_slider = Slider(
    fig[1, 2],
    range = 1:1:k,
    startvalue = 1
)

Label(fig[1, 3], lift(n -> "n = $n", n_slider.value))

k=100
sel_index = 1
contact_coords = random_closed_sphere_curve()
conn = ConstructConnection(contact_coords,k)
trans_dev = [HolomorphicTransform(conn,i) for i in 1:k]
holomorph3,cone_circles = Holomorph3(conn)

contact_curve = [S2_from_angles(c[1],c[2]) for c in contact_coords]
dev_curve = conn.D

xs_contact = [p[1] for c in contact_curve]
ys_contact = [p[2] for c in contact_curve]
zs_contact = [p[3] for c in contact_curve]

xs_dev = [p.re for p in dev_curve]
ys_dev = [p.im for p in dev_curve]
zs_dev = zeros(k)

xs_trans_dev = [p.re for p in trans_dev[sel_index]]
ys_trans_dev = [p.im for p in trans_dev[sel_index]]
zs_trans_dev = zeros(k)

xs_holomorph = [p[1] for p in holomorph3[sel_index]]
ys_holomorph = [p[2] for p in holomorph3[sel_index]]
zs_holomorph = [p[3] for p in holomorph3[sel_index]]

xs_cone_circle = [p[1] for p in cone_circles[sel_index]]
ys_cone_circle = [p[2] for p in cone_circles[sel_index]]
zs_cone_circle = [p[3] for p in cone_circles[sel_index]]

onany(n_slider.value) do
    sel_index = n_slider_value
    xs_trans_dev = [p.re for p in trans_dev[sel_index]]
    ys_trans_dev = [p.im for p in trans_dev[sel_index]]
    zs_trans_dev = zeros(k)
    xs_holomorph = [p[1] for p in holomorph3[sel_index]]
    ys_holomorph = [p[2] for p in holomorph3[sel_index]]
    zs_holomorph = [p[3] for p in holomorph3[sel_index]]
    xs_cone_circle = [p[1] for p in cone_circles[sel_index]]
    ys_cone_circle = [p[2] for p in cone_circles[sel_index]]
    zs_cone_circle = [p[3] for p in cone_circles[sel_index]]
    lines!(ax_trans_dev,xs_trans_dev,ys_trans_dev,zs_trans_dev)
    lines!(ax_holomorph,xs_holomorph,ys_holomorph,zs_holomorph)
    lines!(ax_holomorph,xs_cone_circle,ys_cone_circle,zs_cone_circle)
end
lines!(ax_contact, xs_contact, ys_contact, zs_contact)
lines!(ax_dev,xs_dev,ys_dev,zs_dev)
scatter!(
    ax_d,
    lift(n -> [xs_d[n]], n_slider.value),
    lift(n -> [ys_d[n]], n_slider.value),
    lift(n -> [zs_d[n]], n_slider.value),
    color = :red,
    markersize = 15
)

display(fig)

