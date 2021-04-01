function prediction(xk,w,dt, numDiodes)
    q = xk[1:4]; # Quaternion portion
    b = xk[5:7]; # Bias portion
    c = xk[8:(7+numDiodes)]
    α = xk[(8+numDiodes):(7 + 2*numDiodes)]
    ϵ = xk[(8 + 2*numDiodes):end]

    theta = (norm(w - b)*dt);
    r = (w - b)/norm(w - b);

    qp = qmult(q, [r*sin(theta/2); cos(theta/2)]); 
    
    xn = [qp; b; c; α; ϵ]; # x at next step
    # xn = [qp; b];
   
    alpha = (w - b);
    skew = -1*[0 -alpha[3] alpha[2]; alpha[3] 0 -alpha[1]; -alpha[2] alpha[1] 0];
    R = (I(3) + skew*sin(dt) + (skew^2)*(1 - cos(dt))); # Rodrigues...?

    # A is jacobian of f(x)
    A = [R -dt*I(3); zeros(3,3) I(3)];

    A = [A                           zeros(6, 3 * numDiodes); 
         zeros(3 * numDiodes, 6)     I(3*numDiodes)];

    # # A needs to be updated with partial with respect to epsilon/alpha 
    # A = [R                          -dt*I(3)                 zeros(3, numDiodes)              zeros(3, numDiodes); 
    #      zeros(3,3)                   I(3)                   zeros(3, numDiodes)              zeros(3, numDiodes); 
    #      zeros(numDiodes,3)     zeros(numDiodes,3)               I(numDiodes)         zeros(numDiodes, numDiodes);
    #      zeros(numDiodes,3)     zeros(numDiodes,3)        zeros(numDiodes, numDiodes)               I(numDiodes)];

    return xn, A

end
