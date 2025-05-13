import numpy as np
l_f = 1
l_s = 1
delta = np.pi/12
r_w = 1
F_FL = 350/4;
F_FR = 350/4;
F_RL = 350/4;
F_RR = 350/4;
T_t = 400;
M_z = 0.05;
u = 1/4*np.array([T_t, T_t, T_t, T_t])+np.random.rand(1,4)*20
W_v = np.array([[1000],[1000]])
W_v = W_v.reshape(1,-1)
u = u.T
rate = 0.0000001


B_v = np.array([[1,1,1,1],[(l_f*np.sin(delta)-l_s*np.cos(delta))/r_w,(l_f*np.sin(delta)+l_s*np.cos(delta))/r_w, -l_s/r_w, l_s/r_w]])
B_w = 1/r_w*np.array([1/F_FL,1/F_FR,1/F_RL,1/F_RR])

E_ref = W_v@(B_v@u-np.array([[T_t],[M_z]]))
F_w = B_w@u

J = E_ref*E_ref+F_w*F_w+np.dot(u.T,u)
dJ = 2*B_v.T@W_v.T@W_v@(B_v@u-np.array([[T_t],[M_z]]))
print(J)
r = 0
while np.dot(dJ.T,dJ) > 0.005 and r < 3000:
    dJ = 2*B_v.T@W_v.T@W_v@(B_v@u-np.array([[T_t],[M_z]]))+2*u*(F_w*F_w)
    u = u-rate*dJ
    r += 1


E_ref = W_v@(B_v@u-np.array([[T_t],[M_z]]))
F_w = B_w@u

J = E_ref*E_ref+F_w*F_w+u.T@u
print(J)
print("Iterations: ",r)
print(u)