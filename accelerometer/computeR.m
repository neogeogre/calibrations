function R = computeR(r, p, y)
% computeR applies roll pitch yaq, original is from l-frame to b-frame

R = [cos(p)*cos(y)                          cos(p)*sin(y)                          -sin(p); 
     sin(r)*sin(p)*cos(y) - cos(r)*sin(y)   sin(r)*sin(p)*sin(y) + cos(r)*cos(y)    sin(r)*cos(p);
     cos(r)*sin(p)*cos(y) + sin(r)*sin(y)   cos(r)*sin(p)*sin(y) - sin(r)*cos(y)    cos(r)*cos(p)];

end