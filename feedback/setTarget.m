function setTarget(a,b,c,d,e)
% This functions sets the initial position of the AUV and the target
% position so that there is no need to initiate the two in the running code
% and the dynamics function

    global p1 p2 q1 q2 v;
    p1=a;p2=b;q1=c;q2=d;v=e;

end

