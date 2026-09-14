const TAngle = Float64
const TC1  = ComplexF64
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
const TPlane3 = Tuple{SVector{3,Float64},SVector{3,Float64}}
const ISO3::TSO3 = @SMatrix[1 0 0;0 1 0;0 0 1]
const ISE3::TSE3 = @SMatrix[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]
const Angles :: Vector{TAngle} = collect(range(-pi,pi,length=360))

const CR = sin(pi/5)
const RA::TAngle = pi/2

Angle(v::Float64) = clamp(v,-pi,pi) :: TAngle
Angle(x::Float64,y::Float64) = angle(Complex(x,y)) :: TAngle
function Plane3(azi::TAngle,polar::TAngle) :: TPlane3
  u = @SVector [-sin(azi),cos(azi),0.0]
  v = @SVector [-(sin(polar) * cos(azi)),-(sin(polar) * sin(azi)),cos(polar)]
  return u, v
end

S1(a::TAngle) = TS1(cos(a),sin(a))
S1(r::Float64,a::TAngle) = TS1(r*cos(a),r*sin(a))

function S2(azi::TAngle,polar::TAngle) :: TS2
    pr=sin(polar)
    x=pr*cos(azi)
    y=pr*sin(azi)
    z=cos(polar)
    return TS2(x,y,z)
end

function S2(r::Float64,azi::Float64,polar::Float64) :: TS2
  pr=sin(polar)
  x=r*pr*cos(azi)
  y=r*pr*sin(azi)
  z=r*cos(polar)
  return TS2(x,y,z)
end

function S3(azi::Float64,polar::Float64,orb::Float64) :: TC2
  z1=Complex(cos(orb),sin(orb))*cos(azi/2)
  z2=Complex(cos(orb-polar),sin(orb-polar))*sin(azi/2)
  return TC2(z1,z2)
end

function S3_R3(p::TC2) :: TR3
  z1,z2 = p[1],p[2]
  d=1-z2.im
  x,y,z=(2*z1.re)/d,(2*z1.im)/d,(2*z2.re)/d
  return TR3(x,y,z)
end

function S3_R3(azi::Float64,polar::Float64,orb::Float64) :: TR3
  z1,z2 = S3(azi,polar,orb)
  return S3_R3(TC2(z1,z2))
end

function StereoProj(x::Float64,y::Float64,z::Float64) :: Tupple{ComplexF64,ComplexF64}
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


function SO3(axis::SVector{3,Float64},ang::Float64) :: TSO3
  x,y,z=axis[1],axis[2],axis[3]
  c,s=cos(ang),sin(ang)
  d=1-c
  return @SMatrix[(c+x^2*d) (x*y*d)-(z*s) (x*z*d)+(y*s);(y*x*d)+(z*s) c+(y^2*d) (y*z*d)-(x*s)]
end

function SO3(azi::Float64,polar::Float64,ang::Float64) :: TSO2
  axis=S2(azi,polar)
  return SO3(axis,ang)
end

function MöbiusRotCoef(v::ComplexF64) :: Tupple{ComplexF64,ComplexF64,ComplexF64,ComplexF64}
  m =abs2(v)
  dn=sqrt(1+m)
  a,b,c=1/dn,v/dn,-conj(v)/dn
  d=a
  return a,b,c,d
end

function HolomorphicTransform(dev::Vector{ComplexF64},n::Int) :: Vector{ComplexF64}
  ref = dev[n]
  a,b,c,d = MöbiusRotCoef(ref)
  return [((a*z)+b)/((c*z)+d) for z in conn.D]
end

function ConeCircle(azi::Float64,polar::Float64,angles::Vector{TAngle}=Angles) :: Vector{TS2}
  u,v = Plane3(azi,polar)
  c=0.8*S2(azi,polar)
  return [c+CR*cos(t)*u+CR*sin(t)*v for t in angles]
end

function ConeCircle(p::TS2,angles::Vector{TAngle}) :: Vector{TS2}
 azi=Angle(p[1],p[2])
 pol=acos(p[3])
 u,v = Plane3(azi,pol)
 c = 0.8*p
 return [c+CR*cos(t)*u+CR*sin(t)*v for t in angles]
end

function TorsionAngle(T1::SVector{3,Float64},T2::SVector{3,Float64},T3::SVector{3,Float64},T4::SVector{3,Float64}) :: Float64
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


function HopfFibre(azi::Float64,polar::Float64,res::Int=360) :: THopfFibre
  ts=res==360 ? Angles : range(-pi,pi,length=res)
  return THopfFibre([TC2(az,polar,t) for t in ts])
end

HopfFibre(v::TS2,res::Int=360) = HopfFibre(Angle(v[1],v[2]),acos(v[3]),res) :: THopfFibre
HopfLink(a::TS2,b::TS2,res::Int=360) = THopfLink(HopfFibre(a,res),HopfFibre(b,res)) :: THopfLink
HopfLink(a::TS2,tor::Float64,res::Int=360) = THopfLink(HopfFibre(a,res),HopfFibre(SO3(a,tor)*a,res)) :: THopfLink


function SE3(m::TSO3,p::SVector{3,Float64}) :: TSE3
  return TSE3([m[1,1] m[1,2] m[1,3] p[1]; m[2,1] m[2,2] m[2,3] p[2]; m[3,1] m[3,2] m[3,3] p[3]; 0 0 0 1])
end

function ΔO(C1::SVector{2,Float64},C2::SVector{2,Float64}) :: TSO3
  VC1 = S2(C1[1],C1[2])
  VC2 = S2(C2[1],C2[2])
  h=normalize(cross(VC1,VC2))
  psi=acos(dot(VC1,VC2))
  return SO3(normalize(cross(p[1],p[2])),acos(dot(p[1],p[2])))
end

function ΔP(O::TSO3,psi::Float64) :: TR3
  vec = TR3([O[3,1],O[3,2],O[3,3]])
  return (psi / normalize(vec)) * vec
end


function adjrange(k::Int,li::Int)::Vector{Tuple{Int,Int}}
  return [(mod1(n,li),mod1(n+1,li)) for n in 1:k]
end

function adjpairs(coll::Vector{TS2}, k::Int) :: Vector{Tuple{TS2,TS2}}
  return map(ip-> (coll[ip[1]],coll[ip[2]]),adjrange(k,lastindex(coll)))
end

function adjpairs(coll::Vector{TSO3},k::Int) :: Vector{Tuple{TSO3,TSO3}}
  return map(ip-> (coll[ip[1]],coll[ip[2]]),adjrange(k,lastindex(coll)))
end

function adjpairs(coll::Vector{TR3}, k::Int) :: Vector{Tuple{TR3,TR3}}
  return map(ip-> (coll[ip[1]],coll[ip[2]]),adjrange(k,lastindex(coll)))
end

function adjpairs(coll::Vector{TC2}, k::Int) :: Vector{Tuple{TC2,TC2}}
  return map(ip-> (coll[ip[1]],coll[ip[2]]),adjrange(k,lastindex(coll)))
end

function adjpairs(coll::Vector{TC1}, k::Int) :: Vector{Tuple{TC1,TC1}}
  return map(ip-> (coll[ip[1]],coll[ip[2]]),adjrange(k,lastindex(coll)))
end

function adjpairs(coll::Vector{Float64}, k::Int) :: Vector{Tuple{Float64,Float64}}
  return map(ip-> (coll[ip[1]],coll[ip[2]]),adjrange(k,lastindex(coll)))
end



function ΔO(c2::SVector{2,Float64},c1::SVector{2,Float64}) :: TSO3
  Vc1 = S2(c1[1],c1[2])
  Vc2 = S2(c2[1],c2[2])
  h=normalize(cross(Vc1,Vc2))
  psi=acos(dot(Vc1,Vc2))
  return SO3(h,psi)
end

function fourier_series(a::Vector{Float64},b::Vector{Float64}) :: Vector{Float64}
  li = lastindex(a)
  tau = 2pi
  if li == 0
    throw(ErrorException("coef lists cannot be empty"))
  end
  if li != lastindex(b)
    throw(ErrorException("coef lists must be of equal size"))
  end
  #res = lastindex(angles)
  out = zeros(Float64,360)
  for h in 1:li
    out .+= (a[h] .* cos.(Angles .* h) .+ (b[h] .* sin.(Angles .* h)))
  end
  return out
end

function sphere_curve(a1::Vector{Float64},a2::Vector{Float64},b1::Vector{Float64},b2::Vector{Float64}) :: Vector{Tuple{TAngle,TAngle}}
 harms = lastindex(a1)
 if harms == 0
   throw(ErrorException("coef list cannot be empty"))
 end
 if !allequal(lastindex,[a1,a2,b1,b2])
   throw(ErrorException("coef lists cannot be empty"))
 end
 azi = fourier_series(a1,a2)
 pol = fourier_series(b1,b2)
 #print(azi)
 #print(pol)
 return [(azi[n],pol[n]) for n in 1:360]
end

CInv(v::ComplexF64) = abs(v) <= 1.0 ? v : 1 / conj(v)






