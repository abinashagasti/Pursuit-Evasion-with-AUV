function y = coordinate_transform(x)
% Coordinate transformation of the state space to linearized space
    c=0.02567;
    y=[x(1);x(4)*cos(x(3))-x(5)*sin(x(3));x(3);c*x(6);x(2);x(5)];
end

