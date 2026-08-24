using StaticArrays
using LinearAlgebra
using Random
using GLMakie
include("Constructors.jl")

disc_pt(z::ComplexF64) = abs(z) <= 1 ? z : 1/conj(z)
circ = [[cos(t),sin(t),0] for t in range(0,2pi,length=360)]
curve_pts = [[rand(-5.0:5.0),rand(-5.0:5.0),0.0] for _ in 1:20]
curve = bezier_curve(curve_pts,360)
cmplx_curve = map(p-> disc_pt(p[1]+p[2]im),curve)
scale = tan(acos(4/5)/2)
fig = Figure(size = (800, 800))
ax1 = Axis3(fig[1,1])
ax2 = Axis3(fig[1,2])
sel_n = 1

function lift2(p::ComplexF64)
  m = abs2(p)
  return [(2*p.re)/(1+m),(2*p.im)/(1+m),(1-m)/(1+m)]
end

function htrans(pts::Vector{ComplexF64},refp::ComplexF64)::Vector{ComplexF64}
  m = abs2(refp)
  dnum = sqrt(1+m)
  a,b,c,d = 1/dnum,refp/dnum,-conj(refp)/dnum,1/dnum
  return map(z-> ((a*z)+b)/((c*z)+d),pts)
end

Dn = curve[sel_n][1]+curve[sel_n][2]im
S = lift2(Dn)
circ = ConeCircle(acos(S[3]),atan(S[2]/S[1]))
scaled_curve = map(z-> scale*z,cmplx_curve)
trans_curve = htrans(scaled_curve,Dn)
lift_trans_curve = map(z-> lift2(z),trans_curve)


zs0 = zeros(Float64,360)
xs1 = [p[1] for p in curve]
ys1 = [p[2] for p in curve]
xs2 = Observable([p[1] for p in circ])
ys2 = Observable([p[2] for p in circ])
zs2 = Observable([p[3] for p in circ])
xs3 = Observable([p[1] for p in lift_trans_curve])
ys3 = Observable([p[2] for p in lift_trans_curve])
zs3 = Observable([p[3] for p in lift_trans_curve])


n_slider = Slider(fig[2,1:2],range = 1:360,startvalue = 1)
Label(fig[2,2:3],lift(n -> "n = $n",n_slider.value),fontsize = 18)
lines!(ax1,xs1,ys1,zs0)
lines!(ax2,xs2,ys2,zs2)
lines!(ax2,xs3,ys3,zs3)
scatter!(ax1,lift(n -> Point3f(xs1[n],ys1[n],zs0[n]),n_slider.value),color = :red,markersize = 15)
scatter!(ax2,lift(n -> Point3f(S[1],S[2],S[3]),n_slider.value),color = :red,markersize = 15)


on(n_slider.value) do n
 sel_n = n
 Dn = curve[n][1]+curve[n][2]im
 S = lift2(Dn)
 circ = ConeCircle(acos(S[3]),atan(S[2]/S[1]))
 trans_curve = htrans(scaled_curve,Dn)
 lift_trans_curve = map(z-> lift2(z),trans_curve)
 xs2 = [p[1] for p in circ]
 ys2 = [p[2] for p in circ]
 zs2 = [p[3] for p in circ]
 xs3 = [p[1] for p in lift_trans_curve]
 ys3 = [p[2] for p in lift_trans_curve]
 zs3 = [p[3] for p in lift_trans_curve]
end

display(fig)
