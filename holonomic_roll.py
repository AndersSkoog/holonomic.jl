import numpy as np
from constants import MIN_VAL

def s2_to_r3(theta,phi):
  pr=np.cos(phi)
  x,y,z = pr*np.cos(theta),pr*np.sin(theta),np.sin(phi)
  return np.array([x,y,z])

def z_from_orient(o:(float,float,float,float),R=1.0):
  v = np.array([0.0,0.0,1.0])
  rot_v = quat_rotate(o,v)
  return R * rot_v[3] if R != 1.0 else rot_v[3]


def random_closed_sphere_curve(n=360, k=5, r=0.1, seed=None):
    rng = np.random.default_rng(seed)
    t = np.linspace(0, 2*np.pi, n)
    theta = np.zeros(n)
    phi = np.zeros(n)

    for i in range(1, k+1):
        A_theta, B_theta = rng.normal(size=2)
        A_phi, B_phi = rng.normal(size=2)
        phase_theta = rng.uniform(0, 2*np.pi)
        phase_phi = rng.uniform(0, 2*np.pi)

        theta += A_theta * np.cos(i*t + phase_theta) + B_theta * np.sin(i*t + phase_theta)
        phi   += A_phi   * np.cos(i*t + phase_phi)   + B_phi   * np.sin(i*t + phase_phi)

    return theta, phi


def rolltranslation(sphere_curve,index,contact,R):
  r=1.0
  #index_next = index + 1 if index < len(sphere_curve[0])-1 else 0
  p1 = s2_to_r3(sphere_curve[0][index],sphere_curve[1][index])
  p2 = s2_to_r3(sphere_curve[0][index+1],sphere_curve[1][index+1])
  ang = np.arccos(np.clip(np.dot(p1,p2), -1.0, 1.0))
  rot_axis = np.cross(p1,p2) / np.linalg.norm(np.cross(p1,p2))
  x,y,z = rot_axis
  c,s,C = np.cos(ang),np.sin(ang),1-np.cos(ang)
  R_inc = np.array([
        [c+x*x*C,x*y*C-z*s,x*z*C+y*s],
        [y*x*C+z*s,c+y*y*C,y*z*C-x*s],
        [z*x*C-y*s,z*y*C+x*s,c+z*z*C]
    ])
  R_new = R @ R_inc
  u,v,n = R_new[:, 0],R_new[:, 1],R_new[:, 2]
  d = np.array([0.0, 0.0, -1.0])
  move_dir = np.cross(n, d)
  move_dir_norm = np.linalg.norm(move_dir)
  if move_dir_norm < MIN_VAL:
    disp = np.zeros(3)
  else:
    move_dir = move_dir / move_dir_norm
    disp = (r * ang) * move_dir

  contact_new = contact + disp
  return contact_new,R_new


def inital_R(sphere_curve):
 sc=sphere_curve
 p1,p2 = s2_to_r3(sc[0][0],sc[0][1]),s2_to_r3(sc[1][0],sc[1][1])
 n1,n2 = p1 / np.linalg.norm(p1), p2 / np.linalg.norm(p2)
 t = n2 - n1
 u = t - np.dot(t,n1) * n1
 v = np.cross(n1,u)
 return np.column_stack((u,v,n1))


"""
if __name__ == "__main__":
  import matplotlib.pyplot as plt
  from mpl_toolkits.mplot3d import Axes3D

  #theta,phi = random_closed_sphere_curve(n=n,seed=16,k=4) # sphere curve
  #S = [np.array([np.sin(phi[i])*np.cos(theta[i]),np.sin(phi[i])*np.sin(theta[i]),np.cos(phi[i])]) for i in range(n)]
  #p = np.array([0.0,0.0,0.0])
  #o = (1.0,0.0,0.0,0.0)

  def construct_disc_curve(theta,phi,r):
    last_o = (1.0,0.0,0.0,0.0)
    last_p = np.array([0.0,0.0,0.0])
    n = len(theta)
    pts = []
    orients = []
    for i in range(1,n):
      p,o = rolltranslation(theta,phi,i,last_o,last_p,r)
      #print(p)
      z = complex(p[0],p[1])
      pts.append(boundaryless_disc(z))
      orients.append(o)
      last_p = p
      last_o = o
    return pts,orients


  #for j in range(len(CD)):
  #  z = z_from_orient(O[j])
  #  x,y = D[j][0],D[j][1]
  #  M.append(np.array([x,y,z]))


  n,r = 360,0.1
  ang = np.linspace(0,tau,n)
  theta,phi = random_closed_sphere_curve(n=n,seed=16,k=4,r=r) # sphere curve
  circle = [np.array([cos(a),sin(a)]) for a in ang]

  #print(theta[0],phi[0])
  #print(theta[-1],phi[0])
  #sx = r * np.sin(phi) * np.cos(theta)
  #sy = r * np.sin(phi) * np.sin(theta)
  #sz = r * np.cos(phi)

  pts,orients = construct_disc_curve(theta,phi,r)
  #pts2 = [np.array([p[0],-p[1],0.0]) for p in non_dup_reverse_array(pts)]
  #pts3 = [np.array([-p[0],p[1],0.0]) for p in non_dup_reverse_array(pts)]
  #pts4 = [np.array([-p[0],-p[1],0.0]) for p in non_dup_reverse_array(pts)]
  #disc_curve = make_closed(pts)
  xs1,ys1,zs1 = zip(*pts)
  cx,cy = zip(*circle)
  #xs2,ys2,zs2 = zip(*pts2)
  #xs3,ys3,zs3 = zip(*pts3)
  #xs4,ys4,zs4 = zip(*pts4)
  fig,ax1 = plt.subplots()
  #ax1 = fig.add_subplot(121, projection='2d')
  ax1.plot(xs1, ys1,'b-',linewidth=0.3)
  ax1.plot(cx,cy,'b-',linewidth=0.3)
  #ax1.plot(xs2, ys2,'b-',linewidth=1)
  #ax1.plot(xs3, ys3,'b-',linewidth=1)
  #ax1.plot(xs4, ys4,'b-',linewidth=1)
  ax1.set_aspect('equal')
  #ax2 = fig.add_subplot(122,projection="3d")
  #ax2.plot(sx,sy,sz,'-b',linewidth=1)
  plt.show()

"""













































  


  
