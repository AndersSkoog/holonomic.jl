using StaticArrays
using LinearAlgebra
using Random
using GLMakie
using Purses
include("lib.jl")
coef_a1 = [0.72,0.12,0.47,0.92]
coef_a2 = [0.24,0.25,0.05,0.24]
coef_b1 = [0.11,0.83,0.33,0.55]
coef_b2 = [0.82,0.12,0.53,0.82]
angles_360  = Purse(collect(range(-pi,pi,360)))
contact_curve = Purse(sphere_curve(coef_a1,coef_a2,coef_b1,coef_b2))
contact_curve_adj = Purse(adjpairs(contact_curve[1:360],360))
contact_pts = Purse([S2(v[1],v[2]) for v in contact_curve[1:360]])
contact_pts_adj = adjpairs(contact_pts[1:360],360)
psi = Purse([acos(dot(p[1],p[2])) for p in contact_pts_adj[1:360]])
conn_delta_orient = Purse([ΔO(cp[1],cp[2]) for cp in contact_curve_adj[1:360]])
conn_orient = Purse(accumulate(*,conn_delta_orient[1:360],init=ISO3))
conn_delta_pos = Purse([ΔP(conn_orient[n],psi[n]) for n in 1:360]))
conn_pos = Purse(accumulate(+,conn_delta_pos[1:360],init=@SVector [0.0,0.0,0.0]))
conn_dev = Purse([Complex(p[1],p[2]) for p in conn_pos[1:360]])
conn_tor = Purse([[conn_pos[n][1],conn_pos[n][2],(conn_orient[n] * @SVector[0.0,0.0,1.0])[3]] for n in 1:360])
conn_torang = Purse([TorsionAngle(conn_tor[n][1],conn_tor[n+1],conn_tor[n+2],conn_tor[n+3]) for n in 1:357])
refS = Purse(InvStereoProj.(conn_dev[1:360]))
circ = Purse(ConeCircle.(refS[1:360],angles_360[1:360]))
print(lastindex(psi))
print(lastindex(conn_delta_orient))
print(lastindex(conn_orient))
print(lastindex(conn_delta_pos))
print(lastindex(conn_pos))
print(lastindex(conn_dev))
print(lastindex(conn_tor))
print(lastindex(conn_torang))
print(lastindex(refS))
print(lastindex(circ))
