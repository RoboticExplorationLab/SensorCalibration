import numpy as np 
# import ulab.numpy as np 

def q_from_DCM(R):
    T = R[0,0] + R[1,1] + R[2,2]
    if (T > R[0,0]) and (T > R[1,1]) and (T>R[2,2]):
        q4 = .5*np.sqrt(1+T)
        r  = .25/q4
        q1 = (R[2,1] - R[1,2])*r
        q2 = (R[0,2] - R[2,0])*r
        q3 = (R[1,0] - R[0,1])*r
    elif R[0,0]>R[1,1] and R[0,0]>R[2,2]:
        q1 = .5*np.sqrt(1-T + 2*R[0,0])
        r  = .25/q1
        q4 = (R[2,1] - R[1,2])*r
        q2 = (R[0,1] + R[1,0])*r
        q3 = (R[0,2] + R[2,0])*r
    elif R[1,1]>R[2,2]:
        q2 = .5*np.sqrt(1-T + 2*R[1,1])
        r  = .25/q2
        q4 = (R[0,2] - R[2,0])*r
        q1 = (R[0,1] + R[1,0])*r
        q3 = (R[1,2] + R[2,1])*r
    else:
        q3 = .5*np.sqrt(1-T + 2*R[2,2])
        r  = .25/q3
        q4 = (R[1,0] - R[0,1])*r
        q1 = (R[0,2] + R[2,0])*r
        q2 = (R[1,2] + R[2,1])*r
    
    q = np.array([q1, q2, q3, q4]);    

    if q4<0:
        q = -q
    
    return q

def dcm_from_q(quat):    # How do i call "hat" here?
    # Takes in a quaternion (scalar last - [q⃗, q]) and returns the DCM  (From Kevin)
    
    v = quat[0:3]
    s = quat[3]

    DCM = np.eye(3) +  2*hat(v) @ (s*np.eye(3) + hat(v))

    return DCM

def hat(v):
    vhat = np.array([[0,   -v[2],  v[1]],
                        [v[2],   0,  -v[0]],
                        [-v[1],  v[0],  0]])

    return vhat; 

def qmult(q, p):
    q0, p0 = q[3], p[3]
    q⃗ = q[0:3]
    p⃗ = p[0:3]

    r0 = q0 * p0 - np.dot(q⃗, p⃗)
    r⃗ = (q0 * p⃗) + (p0 * q⃗) + np.cross(q⃗, p⃗)

    r = np.hstack([r⃗, r0])
    return r

def qconj(q):
    q0 = q[3]
    q⃗ = q[0:3]
    return np.hstack([-q⃗, q0])
