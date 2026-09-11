
coef_a1 = [0.72,0.12,0.47,0.92]
coef_a2 = [0.24,0.25,0.05,0.24]
coef_b1 = [0.11,0.83,0.33,0.55]
coef_b2 = [0.82,0.12,0.53,0.82]
angles_360  = Purse(collect(range(-pi,pi,360)))
contact_curve = Purse(sphere_curve(coef_a1,coef_a2,coef_b1,coef_b2,angles_360[1]))
contact_pts = Purse([S2(v[1],v[2]) for v in contact_curve[1]])
contact_pts_adj = Purse(adjpairs(contact_pts[1],360))
psi = Purse([acos(dot(p[1],p[2])) for p in contact_pts_adj])
conn_delta_orient = Purse([ΔO(cp[1],cp[2])) for cp in contact_curve_adj])
conn_orient = Purse(accumulate(*,conn_delta_orient[1],init=ISO3))
conn_delta_pos = Purse([ΔP(O[1][n],psi[1][n]) for n in 1:360]))
conn_pos = Purse(accumulate(+,conn_delta_pos[1],init=@SVector [0.0,0.0,0.0]))
conn_dev = Purse([Complex(p[1],p[2]) for p in conn_pos[1]])
conn_tor = Purse([[conn_pos[1][n][1],conn_pos[1][n][1],(conn_orient[1][n] * @SVector[0.0,0.0,1.0])[3]] for n in 1:360])
conn_torang = Purse([TorsionAngle(conn_tor[1][n],conn_tor[1][n+1],conn_tor[1][n+2],conn_tor[1][n+3]) for n in 1:357])
refS = Purse(InvStereoProj.(conn_dev[1]))
circ = Purse(ConeCircle.(refS[1],angles_360[1]))
holom2 = Purse([HolomorphicTransform(conn_dev[1][n],n) for n in 1:360])
holom3 = Purse([-InvStereoProj.(curve) for curve in holom2[1]]))
#holom4_V1 = Purse([HopfFibre.(curve) for curve in holom3[1]])
#holom4_V2 = Purse([HopfLink.(holom3[1][n],conn_torang[1][n]) for n in 1:357)])
#proj_holom4_V1 = Purse([S3_R3.(fibre) for fibre in holom4_V1[1]])
#proj_holom4_V2 = Purse([S3_R3.(fibre) for fibre in holom4_V1[1]])


