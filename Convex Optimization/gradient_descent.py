import numpy as np
import time

l_f = 1
l_s = 1
delta = np.pi/12
r_w = 1
F_FL = 350/4;
F_FR = 350/4;
F_RL = 350/4;
F_RR = 350/4;
motorLimit = np.array([15*13.5,15*13.5,130*3.4])
t = 8;
T_t = 647;
M_z = 3;
u = 1/4*np.array([T_t, T_t, T_t, T_t]).reshape(1,-1);
W_v = np.array([[1000],[1000]])
W_v = W_v.reshape(1,-1)
u = u.T
rate = 0.00000001

# start = time.time()
B_v = np.array([[1,1,1,1],[(l_f*np.sin(delta)-l_s*np.cos(delta))/r_w,(l_f*np.sin(delta)+l_s*np.cos(delta))/r_w, -l_s/r_w, l_s/r_w]])
B_w = 1/r_w*np.array([1/F_FL,1/F_FR,1/F_RL,1/F_RR]).reshape(1,-1)
c = np.array([[T_t],[M_z]]).reshape(-1,1)

E_ref = W_v@(B_v@u-c)
F_w = B_w@u


J = E_ref*E_ref+F_w*F_w+np.dot(u.T,u)
dbarrier = np.array([1/(-u[0]+motorLimit[0])-1/(u[0]+motorLimit[0]),1/(-u[1]+motorLimit[1])-1/(u[1]+motorLimit[1]),1/(u[2]+u[3]-motorLimit[2]),1/(u[2]+u[3]-motorLimit[2])])
dJ = 2*(B_v.T@W_v.T@W_v@(B_v@u-np.array([[T_t],[M_z]]))+B_w.T*np.dot(B_w,u)+u)+0*dbarrier;
print(J)
# A = u.T@(B_v.T@W_v.T@W_v@B_v+B_w.T@B_w+np.identity(4))
# B =B_v.T@W_v.T@W_v@c
# B = B.T
r = 0
print(dJ)
# while np.dot(dJ.T,dJ) > 0.005 and r < 30000:
#     dbarrier = np.array([1/(-u[0]+motorLimit[0])-1/(u[0]+motorLimit[0]),1/(-u[1]+motorLimit[1])-1/(u[1]+motorLimit[1]),1/(u[2]+u[3]-motorLimit[2]),1/(u[2]+u[3]-motorLimit[2])])
#     dJ = 2*B_v.T@W_v.T@W_v@(B_v@u-np.array([[T_t],[M_z]]))+2*B_w.T*np.dot(B_w,u)+u+dbarrier
#     u = u-rate*dJ
#     r += 1

# end = time.time()
# E_ref = W_v@(B_v@u-np.array([[T_t],[M_z]]))
# F_w = B_w@u

# J = E_ref*E_ref+F_w*F_w+u.T@u
# print(J)
# print("Iterations: ",r)
# print(u)
# print("Time:",end-start)