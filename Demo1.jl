using StaticArrays
using LinearAlgebra
using Random
using GLMakie
include("Constructors.jl")

k = 100

contact_coords =
    random_closed_sphere_curve()

conn =
    ConstructConnection(
        contact_coords,
        k
    )

# ------------------------------------------------------------
# Contact curve
# ------------------------------------------------------------

contact_curve = S2[
    S2_from_angles(
        c[1],
        c[2]
    )
    for c in contact_coords
]

xs_contact,
ys_contact,
zs_contact = xyz(contact_curve)

# ------------------------------------------------------------
# Development
# ------------------------------------------------------------

dev_curve = conn.D

xs_dev = [
    z.re
    for z in dev_curve
]

ys_dev = [
    z.im
    for z in dev_curve
]

zs_dev = zeros(Float64, length(dev_curve))

# ============================================================
# Initial selected transformation
# ============================================================

n0 = 1

trans0 =
    HolomorphicTransform(
        conn,
        n0
    )

H30,
circle0 =
    Holomorph3(
        conn,
        n0
    )

# ------------------------------------------------------------
# Dynamic transformed-development data
# ------------------------------------------------------------

trans_x = Observable([
    z.re
    for z in trans0
])

trans_y = Observable([
    z.im
    for z in trans0
])

trans_z = Observable(
    zeros(Float64, length(trans0))
)

# ------------------------------------------------------------
# Dynamic lifted curve
# ------------------------------------------------------------

holo_x = Observable([
    p[1]
    for p in H30
])

holo_y = Observable([
    p[2]
    for p in H30
])

holo_z = Observable([
    p[3]
    for p in H30
])

# ------------------------------------------------------------
# Dynamic cone circle
# ------------------------------------------------------------

circle_x = Observable([
    p[1]
    for p in circle0
])

circle_y = Observable([
    p[2]
    for p in circle0
])

circle_z = Observable([
    p[3]
    for p in circle0
])

# ============================================================
# Figure
# ============================================================

fig = Figure(
    size = (1600, 850)
)

# ------------------------------------------------------------
# Titles
# ------------------------------------------------------------

Label(
    fig[1, 1],
    "Contact curve",
    fontsize = 20
)

Label(
    fig[1, 2],
    "Development",
    fontsize = 20
)

Label(
    fig[1, 3],
    "Transformed development",
    fontsize = 20
)

Label(
    fig[1, 4],
    "Lifted to sphere",
    fontsize = 20
)

# ------------------------------------------------------------
# Axes
# ------------------------------------------------------------

ax_contact =
    Axis3(fig[2, 1])

ax_dev =
    Axis3(fig[2, 2])

ax_trans_dev =
    Axis3(fig[2, 3])

ax_holomorph =
    Axis3(fig[2, 4])

# ============================================================
# Static plots
# ============================================================

lines!(
    ax_contact,
    xs_contact,
    ys_contact,
    zs_contact,
    linewidth = 2
)

lines!(
    ax_dev,
    xs_dev,
    ys_dev,
    zs_dev,
    linewidth = 2
)

# ============================================================
# Dynamic plots
# ============================================================

lines!(
    ax_trans_dev,
    trans_x,
    trans_y,
    trans_z,
    linewidth = 2
)

lines!(
    ax_holomorph,
    holo_x,
    holo_y,
    holo_z,
    linewidth = 2
)

lines!(
    ax_holomorph,
    circle_x,
    circle_y,
    circle_z,
    linewidth = 1
)

# ============================================================
# Slider
# ============================================================

n_slider = Slider(
    fig[4, 2:3],
    range = 1:k,
    startvalue = 1
)

Label(
    fig[3, 2:3],
    lift(
        n -> "n = $n",
        n_slider.value
    ),
    fontsize = 18
)

# ============================================================
# Selected point
# ============================================================

scatter!(
    ax_dev,

    lift(
        n -> Point3f(
            xs_dev[n],
            ys_dev[n],
            zs_dev[n]
        ),
        n_slider.value
    ),

    color = :red,
    markersize = 15
)



# ============================================================
# Slider callback
# ============================================================

on(n_slider.value) do n

    # --------------------------------------------------------
    # Holomorphic transform
    # --------------------------------------------------------

    trans =
        HolomorphicTransform(
            conn,
            n
        )

    trans_x[] = [
        z.re
        for z in trans
    ]

    trans_y[] = [
        z.im
        for z in trans
    ]

    # --------------------------------------------------------
    # Lift to sphere
    # --------------------------------------------------------

    H3, circle =
        Holomorph3(
            conn,
            n
        )

    holo_x[] = [
        p[1]
        for p in H3
    ]

    holo_y[] = [
        p[2]
        for p in H3
    ]

    holo_z[] = [
        p[3]
        for p in H3
    ]

    # --------------------------------------------------------
    # Cone circle
    # --------------------------------------------------------

    circle_x[] = [
        p[1]
        for p in circle
    ]

    circle_y[] = [
        p[2]
        for p in circle
    ]

    circle_z[] = [
        p[3]
        for p in circle
    ]
end

display(fig)
