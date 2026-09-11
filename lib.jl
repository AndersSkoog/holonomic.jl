const TS1  = SVector{2, Float64}
const TS2  = SVector{3, Float64}
const TR3  = SVector{3, Float64}
const TR2  = SVector{2, Float64}
const TC2  = SVector{2, ComplexF64}
const TSO3 = SMatrix{3, 3, Float64, 9}
const TSU2 = SMatrix{2, 2, ComplexF64, 4}
const TSE3 = SMatrix{4, 4, Float64, 16}
const THopfFibre = Vector{TC2}
const THopfLink = Tuple{THopfFibre,THopfFibre}
const ISO3::TSO3 = @SMatrix [1 0 0;0 1 0;0 0 1]
const ISE3::TSE3 = @SMatrix [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]
const default_angles = range(-pi,pi,length=360)
const CR = sin(pi/5)
const RA = pi/2

Angle(v::Float64,direction::Int) = [clamp(v,-pi,pi),clamp(v,0,2pi)][direction]
Angle(x::Float64,y::Float64) = angle(Complex(x,y))

function Plane3(azi::Float64,polar::Float64)
  u = @SVector [-sin(azi),cos(azi),0]
  v = @SVector [-(sin(polar) * cos(azi)),-(sin(polar) * sin(azi)),cos(polar))
  return u, v
end

S1(a::Float64) = TS1(cos(a),sin(a))
S1(r::Float64,a::Float64) = TS1(r*cos(a),r*sin(a))

function S2(azi::Float64,polar::Float64)
    pr=sin(polar)
    x=pr*cos(azi)
    y=pr*sin(azi)
    z=cos(polar)
    return TS2(x,y,z)
end

function S2(r::Float64,azi::Float64,polar::Float64)
  pr=sin(polar)
  x=r*pr*cos(azi)
  y=r*pr*sin(azi)
  z=r*cos(polar)
  return TS2(x,y,z)
end

function S3(azi::Float64,polar::Float64,orb::Float64)
  z1=Complex(cos(orb),sin(orb))*cos(azi/2)
  z2=Complex(cos(orb-polar),sin(orb-polar))*sin(azi/2)
  return TC2(z1,z2)
end

function S3_R3(p::TC2)
  z1,z2 = p[1],p[2]
  d=1-z2.im
  x,y,z=(2*z1.re)/d,(2*z1.im)/d,(2*z2.re)/d
  return TR3(x,y,z)
end

function S3_R3(azi::Float65,polar::Float64,orb::Float64)
  z1,z2 = S3(azi,polar,orb)
  return S3_R3(TC2(z1,z2))
end

function StereoProj(x::Float64,y::Float64,z::Float64) :: ComplexF64
  ζ = (x + y * im) / (1 - z)
  ξ = (x - y * im) / (1 + z)
  return ζ,ξ
end

StereoProj(v::TS2) = StereoProj(v[1],v[2],v[3])

function InvStereoProj(v::SVector{2,Float64}) :: TS2
  m = abs2(v)
  x = (2*v.re)/(1+m)
  y = (2*v.im)/(1+m)
  z = (1-m)/(1+m)
  return TS2(x,y,z)
end

function InvStereoProj(v::ComplexF64) :: TS2
    m = abs2(v)
    x = (2*v.re)/(1+m)
    y = (2*v.im)/(1+m)
    z = (1-m)/(1+m)
    return TS2(x,y,z)
end


function SO3(axis::SVector{3,Float64},ang::Float64)
  x,y,z=axis[1],axis[2],axis[3]
  c,s=cos(ang),sin(ang)
  d=1-c
  return @SMatrix [(c+x^2*d) (x*y*d)-(z*s) (x*z*d)+(y*s);(y*x*d)+(z*s) c+(y^2*d) (y*z*d)-(x*s)]
end

function SO3(azi::Float,polar::Float64,ang::Float64)
  axis=S2(azi,polar)
  return SO3(axis,ang)
end

function MöbiusRotCoef(v::ComplexF64)
  m =abs2(v)
  dn=sqrt(1+m)
  a,b,c=1/dn,v/dn,-conj(v)/dn
  d=a
  return a,b,c,d
end

function HolomorphicTransform(dev::Vector{ComplexF64},n::Int)
  ref = dev[n]
  a,b,c,d = MöbiusRotCoef(ref)
  return [((a*z)+b)/((c*z)+d) for z in conn.D]
end

function ConeCircle(azi::Float64,polar::Float64,res::Int = 360)
  ts=res==360 ? default_angles : range(0,2π,length=res)
  u,v = Plane3(azi,polar)
  c=0.8*S2(azi,polar)
  return [@SVector c+CR*cos(t)*u+CR*sin(t)*v for t in ts]
end

function ConeCircle(p::S2,angles:Vector{Float64})
 azi=Angle(p[1],p[2])
 pol=acos(p[3])
 u,v = Plane3(azi,pol)
 c = 0.8*p
 return [@SVector c+CR*cos(t)*u+CR*sin(t)*v for t in angles]
end

function TorsionAngle(T1::SVector{3,Float64},T2::SVector{3,Float64},T3::SVector{3,Float64},T4::SVector{3,Float64})
  d1 = T2-T1
  d2 = T3-T2
  d3 = T4-T3
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
  return Angle(x,y)
end

function TorsionAngle(T::Vector{SVector{3,Float64}},n::Int)::Float64
  return TorsionAngle(T[n],T[n+1],T[n+2],T[n+3])
end


function HopfFibre(azi::Float64,polar::Float64,res::Int=360)
  ts=res==360 ? default_angles : range(-pi,pi,length=res)
  return THopfFibre([TC2(az,polar,t) for t in ts])
end

HopfFibre(v::S2,res::Int=360) = HopfFibre(Angle(v[1],v[2]),acos(v[3]),res)
HopfLink(a::S2,b::S2,res::Int=360) = THopfLink(HopfFibre(a,res),HopfFibre(b,res))
HopfLink(a::S2,tor::Float64,res::Int=360) = THopfLink(HopfFibre(a,res),HopfFibre(SO3(a,tor)*a,res))


function SE3(m::SO3,p::SVector{3,Float64})
  return @SMatrix{4,Float64,16} [m[1,1] m[1,2] m[1,3] p[1]; m[2,1] m[2,2] m[2,3] p[2]; m[3,1] m[3,2] m[3,3] p[3]; 0 0 0 1]
end

function ΔO(C1::SVector{2,Float64},C2::SVector{2,Float64})
  VC1 = S2(C1[1],C1[2])
  VC2 = S2(C2[1],C2[2])
  h=@SVector normalize(cross(VC1,VC2))
  psi=acos(dot(VC1,VC2))
  return SO3(normalize(cross(p[1],p[2])),acos(dot(p[1],p[2])))
end

function ΔP(O::TSO3,psi::Float64)
  vec = @SVector [O[3,1],O[3,2],O[3,3]]
  return (psi / normalize(vec)) * vec
end


function adjrange(k::Int,li::Int)::Vector{Tuple(Int,Int)}
  return [(mod1(n,li),mod1(n+1,li)) for n in 1:k]
end

adjpairs(coll,k) = map(ip-> (coll[ip[1]],coll[ip[2]]),adjrange(k,lastindex(arr)))

function ΔO(C1::SVector{2,Float64},C2::SVector{2,Float64})
  VC1 = S2(C1[1],C1[2])
  VC2 = S2(C2[1],C2[2])
  h=@SVector normalize(cross(VC1,VC2))
  psi=acos(dot(VC1,VC2))
  return SO3(h,psi)
end

function ΔP(O::TSO3,psi::Float64)
  vec = @SVector [O[3,1],O[3,2],O[3,3]]
  return (psi / normalize(vec)) * vec
end

function fourier_series(a::Vector{Float64},b::Vector{Float64},angles::Vector{Float64})
  li = lastindex(cA)
  tau = 2pi
  if li == 0
    throw(ErrorException("coef lists cannot be empty"))
  end
  if li != lastindex(cB)
    throw(ErrorException("coef lists must be of equal size"))
  end
  res = lastindex(angles)
  out = zeroes(Float64,res)
  for h in 1:li
    out .+= (a[h] .* cos.(angles .* h) .+ (b[h] .* sin.(angles .* h)))
  end
  return out
end

function sphere_curve(a1::Vector{Float64},a2::Vector{Float64},b1::Vector{Float64},b2::Vector{Float64},angles::Vector{Float64})
 harms = lastindex(a1)
 if hams == 0
   throw(ErrorException("coef list cannot be empty"))
 end
 if !allequal(lastinex,[a1,a2,b1,b2])
   throw(ErrorException("coef lists cannot be empty"))
 end
 azi = fourier_series(a1,a2,angles)
 pol = fourier_series(b1,b2,angles)
 return [@SVector [v[1],v[2]] v for in zip(azi,pol)]
end

CInv(v::ComplexF64) = abs(v) <= 1.0 ? v : 1 / conj(v)

function xyz(curve::Vector{S2}) :: Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}
  x = [p[1] for p in curve]
  y = [p[2] for p in curve]
  z = [p[3] for p in curve]
  return x, y, z
end






