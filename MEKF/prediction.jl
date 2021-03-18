function prediction(xk,w,dt)
    q = xk[1:4]; # Quaternion portion
    b = xk[5:7]; # Bias portion

    theta = (norm(w - b)*dt);
    r = (w - b)/norm(w - b);

    # qp = qmult(q, [cos(theta/2); r*sin(theta/2)]);  #################
    qp = qmult(q, [r*sin(theta/2); cos(theta/2)]); 
    
    
    xn = [qp; b]; # x at next step
   
    alpha = (w - b);
    skew = -1*[0 -alpha[3] alpha[2]; alpha[3] 0 -alpha[1]; -alpha[2] alpha[1] 0];
    R = (I(3) + skew*sin(dt) + (skew^2)*(1 - cos(dt))); # Rodrigues...?

    A = [R -dt*I(3); zeros(3,3) I(3)];

    return xn, A

end
