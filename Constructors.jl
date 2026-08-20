const S2 = SVector{3,Float64}
const R3 = SVector{3,Float64}
const SO3 = SMatrix{3,3,Float64,9}
const SU2 = SMatrix{2,2,ComplexF64,4}
const su2 = SMatrix{2,2,ComplexF64,4}
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

angle_res = 360
angles = [(2pi/angle_res)*ind for ind in 1:angle_res]
basefiber = [exp(theta*im) for theta in angles]

function S2_from_angles(θ::Float64,φ::Float64)::S2
  x=cos(φ)*cos(θ)
  y=cos(φ)*sin(θ)
  z=sin(φ)
  v=normalize([x,y,z])
  return S2(v[1],v[2],v[3])
end

function SO3_from_axis_angle(axis::S2,ang::Float64)
    x,y,z = axis.x,axis.y,axis.z
    c = cos(ang)
    s = sin(ang)
    C = 1 - c
    R = [
        c + x*x*C   x*y*C - z*s   x*z*C + y*s;
        y*x*C + z*s  c + y*y*C     y*z*C - x*s;
        z*x*C - y*s  z*y*C + x*s   c + z*z*C
        ]
    return SO3(R)
end


#Plane tangent to the unit-sphere coordinate (θ,φ)
function SpherePlane(θ::Float64,φ::Float64)
    c1=sin(θ)
    c2=sin(φ)*cos(θ)
    c3=cos(θ)
    c4=sin(φ)*sin(θ)
    c5=cos(φ)
    u=(-c1,c3,0)
    v=(-c2,-c4,c5)
    return u,v
end

#intersecting-circle of the sphere and the cone with radius 1 at a distance 2 from the origin
function ConeCircle(θ::Float64,φ::Float64,res::Int64=360)
  r=sin(pi/5)
  u,v = SpherePlane(θ,φ)
  center = (4/5)*S2_from_angles(θ,φ)
  return [center + (cos(t) * r * u) + (sin(t) * r * v) for t in range(0,(2*pi)/res,2*pi)]
end

function TorsionAngle(T,n::Int64)
  d1 = T[n+1] - T[n]
  d2 = T[n+2] - T[n+1]
  d3 = T[n+3] - T[n+2]
  n1 = normalize(cross(d1,d2))
  n2 = normalize(cross(d2,d3))
  x = dot(n1,n2)
  y = dot(normalize(d2),cross(n1,n2))
  return atan(y,x)
end

function StereoProj(v::S2)
  x,y,z = v[1],v[2],v[3]
  ζ = (v.x+v.y*im) / (1-v.z)
  ξ = (v.x-v.y*im) / (1+v.z)
  return ζ,ξ
end

function InvStereoProj(v::ComplexF64)
  m = abs2(v)
  x = (2*v.re)/(1+m)
  y = (2*v.im)/(1+m)
  z = (1-m)/(1+m)
  return S2(x,y,z)
end

function random_closed_sphere_curve(n::Int=360, h::Int=5)
  t = range(0, 2π, length=n)
  θ = zeros(n)
  φ = zeros(n)
  for i in 1:h
    Aθ, Bθ = randn(2)
    Aφ, Bφ = randn(2)
    phaseθ = 2π * rand()
    phaseφ = 2π * rand()
    θ .+= Aθ .* cos.(i .* t .+ phaseθ) .+
    Bθ .* sin.(i .* t .+ phaseθ)
    φ .+= Aφ .* cos.(i .* t .+ phaseφ) .+
    Bφ .* sin.(i .* t .+ phaseφ)
  end
  return [(θ[i],φ[i]) for i in 1:n]
end

function HolomorphicTransform(conn::SphereConnection,n::Int64)
  D_n = conn.D[n]
  m=abs2(D_n)
  dnum=sqrt(1+m)
  a=1/dnum
  b=D_n/dnum
  c=-conjugate(D_n)/dnum
  d=a
  return [((a*z)+b)/((c*z)+d) for z in conn.W]
end


function ConstructConnection(contacts::Vector{Tuple{Float64,Float64}},k::Int64) :: SphereConnection
  l = lastindex(contacts)
  down = [0.0,0.0,-1.0]
  s=tan(acos(4/5)/2)
  O_prev = [1 0 0; 0 1 0; 0 0 1]
  Dev_prev = [0,0,0]
  Theta = [contacts[i][1] for i in 1:l]
  Phi = [contacts[i][2] for i in 1:l]
  C = [S2_from_angles(Theta[i],Phi[i]) for i in 1:l]
  Alpha = [cos(Theta[i]/2) + 0im for i in 1:l]
  Beta = [exp(Phi[i]*im)*sin(Theta[i]/2) for i in 1:l]
  O = []
  D = []
  T = []
  S = []
  W = []
  for n in 1:k
    #ind = mod1(n,l)
    #print(mod1(n,l))
    #print(mod1(n+1,l))
    p1 = C[mod1(n,l)]
    p2 = C[mod1(n+1,l)]
    rotaxis = cross(p1,p2)
    ang = acos(dot(p1,p2))
    O_inc = SO3_from_axis_angle(rotaxis,ang)
    O_new = O_prev * O_inc
    normal = [O_new[3,1],O_new[3,2],O_new[3,3]]
    move_dir = cross(normal,down)
    move_norm = norm(move_dir)
    disp = move_norm < 1e-12 ? [0.0,0.0,0.0] : ang * (move_dir / move_norm)
    Dev_new = Dev_prev + disp
    T_vec = O_new * [0.0,0.0,0.1]
    T_n = [Dev_new[1],Dev_new[2],T_vec[3]]
    D_n = Dev_new[1]+Dev_new[2]im
    push!(O,O_new)
    push!(D,D_n)
    push!(T,T_n)
    push!(S,InvStereoProj(D_n))
    push!(W,Beta[n]*(s*D_n))
  end
  A=[TorsionAngle(T,n) for n in 1:k-3]
  R=[SO3_from_axis_angle(S[n],A[n]) for n in 1:k-3]
  return SphereConnection(C,Theta,Phi,Alpha,Beta,O,D,T,S,W,A,R)
end


function HopfLink(conn::SphereConnection,n::Int64)
  a = p
  b = a * R[n]
  theta_a = acos(a[3])
  phi_a = atan(a[2],a[1])
  alpha_a =cos(theta_a/2)+0im
  beta_a = exp(phi_a*im)*sin(theta_a/2)
  theta_b = acos(b[3])
  phi_b = atan(b[2],b[1])
  alpha_b =cos(theta_b/2)+0im
  beta_b = exp(phi_b*im)*sin(theta_b/2)
  z1= [alpha_a beta_a;-conjugate(beta_a) conjugate(alpha_a)] * basefiber
  z2= [alpha_b beta_b;-conjugate(beta_b) conjugate(alpha_b)] * basefiber
  return (z1,z2)
end

function Holomorph3(conn::SphereConnection)
  k = lastindex(conn.D)
  H=[]
  ConeCirc=[]
  for n in 1:k
    m = abs2(D[n])
    dnum = sqrt(1+m)
    a=1/dnum
    b=D[n]/dnum
    c=-conjugate(D[n])/dnum
    d=a
    Circ_n=ConeCircle(acos(conn.S[3]),atan(conn.S[2],conn.S[1]))
    H_n=[InvStereoProj(v) for v in HolomorphicTransform(conn,n)]
    push!(H,H_n)
    push!(ConeCirc,Circ_n)
  end
  return H,ConeCirc
end

function Holomorph4(conn::SphereConnection,n::Int64)
  k = lastindex(conn.D)
  F = []
  for ind in 1:(k-3)
    push!(F,HopfLink(S[n],H[n],TorsionAngle(T[n])))
  end
  return F
end


