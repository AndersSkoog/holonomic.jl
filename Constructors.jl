# ============================================================
# Types
# ============================================================

const S2  = SVector{3, Float64}
const R3  = SVector{3, Float64}
const SO3 = SMatrix{3, 3, Float64, 9}
const SU2 = SMatrix{2, 2, ComplexF64, 4}

struct SphereConnection
    C::Vector{S2}
    θ::Vector{Float64}
    φ::Vector{Float64}
    α::Vector{ComplexF64}
    β::Vector{ComplexF64}
    O::Vector{SO3}
    D::Vector{ComplexF64}
    T::Vector{R3}
    S::Vector{S2}
    W::Vector{ComplexF64}
    A::Vector{Float64}
    R::Vector{SO3}
end


# ============================================================
# Global fiber
# ============================================================

const ANGLE_RES = 360

const angles = range(
    0,
    2π,
    length = ANGLE_RES + 1
)[1:end-1]

const basefiber = ComplexF64[
    exp(im * θ)
    for θ in angles
]


# ============================================================
# S² from spherical angles
# ============================================================

function S2_from_angles(θ::Float64, φ::Float64)::S2

    cφ = cos(φ)

    return S2(
        cφ * cos(θ),
        cφ * sin(θ),
        sin(φ)
    )
end


# ============================================================
# SO(3) rotation from axis-angle
# ============================================================

function SO3_from_axis_angle(
    axis::S2,
    ang::Float64
)::SO3

    x, y, z = axis

    c = cos(ang)
    s = sin(ang)
    C = 1 - c

    return SO3(
        c + x*x*C,
        x*y*C - z*s,
        x*z*C + y*s,

        y*x*C + z*s,
        c + y*y*C,
        y*z*C - x*s,

        z*x*C - y*s,
        z*y*C + x*s,
        c + z*z*C
    )
end


# ============================================================
# Tangent plane to sphere
# ============================================================

function SpherePlane(
    θ::Float64,
    φ::Float64
)

    cθ = cos(θ)
    sθ = sin(θ)
    cφ = cos(φ)
    sφ = sin(φ)

    u = S2(
        -sθ,
         cθ,
         0.0
    )

    v = S2(
        -sφ * cθ,
        -sφ * sθ,
         cφ
    )

    return u, v
end


# ============================================================
# Circle formed by sphere/cone intersection
# ============================================================

function ConeCircle(
    θ::Float64,
    φ::Float64,
    res::Int = 360
)

    r = sin(π / 5)

    u, v = SpherePlane(θ, φ)

    center = (4 / 5) * S2_from_angles(θ, φ)

    return S2[
        center +
        cos(t) * r * u +
        sin(t) * r * v
        for t in range(
            0,
            2π,
            length = res + 1
        )
    ]
end


# ============================================================
# Discrete torsion / dihedral angle
# ============================================================

function TorsionAngle(
    T::Vector{R3},
    n::Int
)::Float64

    d1 = T[n+1] - T[n]
    d2 = T[n+2] - T[n+1]
    d3 = T[n+3] - T[n+2]

    c1 = cross(d1, d2)
    c2 = cross(d2, d3)

    # Degenerate cases
    if norm(c1) < 1e-12 || norm(c2) < 1e-12
        return 0.0
    end

    n1 = normalize(c1)
    n2 = normalize(c2)
    axis = normalize(d2)

    x = dot(n1, n2)
    y = dot(axis, cross(n1, n2))

    return atan(y, x)
end


# ============================================================
# Stereographic projection
# ============================================================

function StereoProj(v::S2)

    x, y, z = v

    ζ = (x + y * im) / (1 - z)
    ξ = (x - y * im) / (1 + z)

    return ζ, ξ
end


# ============================================================
# Inverse stereographic projection
# ============================================================

function InvStereoProj(
    v::ComplexF64
)::S2

    m = abs2(v)

    return S2(
        2 * v.re / (1 + m),
        2 * v.im / (1 + m),
        (1 - m) / (1 + m)
    )
end


# ============================================================
# Random closed sphere curve
# ============================================================

function random_closed_sphere_curve(
    n::Int = 360,
    h::Int = 5
)

    t = range(
        0,
        2π,
        length = n + 1
    )[1:end-1]

    θ = zeros(Float64, n)
    φ = zeros(Float64, n)

    for i in 1:h

        Aθ, Bθ = randn(2)
        Aφ, Bφ = randn(2)

        phaseθ = 2π * rand()
        phaseφ = 2π * rand()

        θ .+=
            Aθ .* cos.(i .* t .+ phaseθ) .+
            Bθ .* sin.(i .* t .+ phaseθ)

        φ .+=
            Aφ .* cos.(i .* t .+ phaseφ) .+
            Bφ .* sin.(i .* t .+ phaseφ)
    end

    return [
        (θ[i], φ[i])
        for i in 1:n
    ]
end


# ============================================================
# Holomorphic transformation
#
# IMPORTANT:
# Only computes the transform for one n.
# ============================================================

function HolomorphicTransform(
    conn::SphereConnection,
    n::Int
)

    D_n = conn.D[n]

    m = abs2(D_n)
    dnum = sqrt(1 + m)

    a = 1 / dnum
    b = D_n / dnum
    c = -conj(D_n) / dnum
    d = a

    return ComplexF64[
        ((a * z) + b) / ((c * z) + d)
        for z in conn.W
    ]
end


# ============================================================
# Construct connection
# ============================================================

function ConstructConnection(
    contacts::Vector{Tuple{Float64, Float64}},
    k::Int
)::SphereConnection

    l = length(contacts)

    down = S2(
        0.0,
        0.0,
       -1.0
    )

    s = tan(acos(4 / 5) / 2)

    O_prev = SO3(I)
    Dev_prev = R3(0.0, 0.0, 0.0)

    # --------------------------------------------------------
    # Input angles
    # --------------------------------------------------------

    θ = [
        contacts[i][1]
        for i in 1:l
    ]

    φ = [
        contacts[i][2]
        for i in 1:l
    ]

    # --------------------------------------------------------
    # Sphere contacts
    # --------------------------------------------------------

    C = S2[
        S2_from_angles(θ[i], φ[i])
        for i in 1:l
    ]

    # --------------------------------------------------------
    # Spinor coordinates
    # --------------------------------------------------------

    α = ComplexF64[
        cos(θ[i] / 2) + 0im
        for i in 1:l
    ]

    β = ComplexF64[
        exp(im * φ[i]) * sin(θ[i] / 2)
        for i in 1:l
    ]

    # --------------------------------------------------------
    # Typed storage
    # --------------------------------------------------------

    O = SO3[]
    D = ComplexF64[]
    T = R3[]
    S = S2[]
    W = ComplexF64[]

    sizehint!(O, k)
    sizehint!(D, k)
    sizehint!(T, k)
    sizehint!(S, k)
    sizehint!(W, k)

    # --------------------------------------------------------
    # Development
    # --------------------------------------------------------

    for n in 1:k

        p1 = C[mod1(n, l)]
        p2 = C[mod1(n + 1, l)]

        d = clamp(dot(p1, p2), -1.0, 1.0)

        ang = acos(d)

        rotaxis = cross(p1, p2)

        # Handle coincident / antipodal points
        if norm(rotaxis) < 1e-12

            O_inc = SO3(I)

        else

            rotaxis = normalize(rotaxis)

            O_inc = SO3_from_axis_angle(
                rotaxis,
                ang
            )
        end

        O_new = O_prev * O_inc

        # ----------------------------------------------------
        # Normal
        # ----------------------------------------------------

        normal = S2(
            O_new[3, 1],
            O_new[3, 2],
            O_new[3, 3]
        )

        # ----------------------------------------------------
        # Development direction
        # ----------------------------------------------------

        move_dir = cross(normal, down)
        move_norm = norm(move_dir)

        disp =
            move_norm < 1e-12 ?
            R3(0.0, 0.0, 0.0) :
            (ang / move_norm) * move_dir

        Dev_new = Dev_prev + disp

        # ----------------------------------------------------
        # Tangent
        # ----------------------------------------------------

        T_vec =
            O_new * @SVector [0.0, 0.0, 0.1]

        T_n = R3(
            Dev_new[1],
            Dev_new[2],
            T_vec[3]
        )

        # ----------------------------------------------------
        # Complex development coordinate
        # ----------------------------------------------------

        D_n =
            Dev_new[1] +
            Dev_new[2] * im

        # ----------------------------------------------------
        # Store
        # ----------------------------------------------------

        push!(O, O_new)
        push!(D, D_n)
        push!(T, T_n)
        push!(S, InvStereoProj(D_n))
        push!(
            W,
            s * D_n
        ) # rotate (s*D_n) to make distinct

        # ----------------------------------------------------
        # IMPORTANT:
        # accumulate the connection
        # ----------------------------------------------------

        O_prev = O_new
        Dev_prev = Dev_new
    end

    # --------------------------------------------------------
    # Torsion
    # --------------------------------------------------------

    A = Float64[
        TorsionAngle(T, n)
        for n in 1:k-3
    ]

    # --------------------------------------------------------
    # Rotations
    # --------------------------------------------------------

    R = SO3[
        SO3_from_axis_angle(
            S[n],
            A[n]
        )
        for n in 1:k-3
    ]

    return SphereConnection(
        C,
        θ,
        φ,
        α,
        β,
        O,
        D,
        T,
        S,
        W,
        A,
        R
    )
end


# ============================================================
# Holomorph 3
#
# Only computes ONE transformed curve and ONE cone circle.
# ============================================================

function Holomorph3(
    conn::SphereConnection,
    n::Int
)

    # --------------------------------------------------------
    # Holomorphic transform
    # --------------------------------------------------------

    H = HolomorphicTransform(
        conn,
        n
    )

    # --------------------------------------------------------
    # Lift transformed complex curve to sphere
    # --------------------------------------------------------

    H3 = S2[
        InvStereoProj(z)
        for z in H
    ]

    # --------------------------------------------------------
    # Cone circle corresponding to S[n]
    # --------------------------------------------------------

    s = conn.S[n]

    θ = acos(
        clamp(s[3], -1.0, 1.0)
    )

    φ = atan(
        s[2],
        s[1]
    )

    circle = ConeCircle(
        θ,
        φ
    )

    return H3, circle
end


# ============================================================
# Convert S2 curve to coordinate arrays
# ============================================================

function xyz(
    curve::Vector{S2}
)

    x = [p[1] for p in curve]
    y = [p[2] for p in curve]
    z = [p[3] for p in curve]

    return x, y, z
end
