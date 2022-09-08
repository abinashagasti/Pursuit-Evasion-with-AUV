function x = reverse_coordinate_transform(y)
% Coordinate transformation of the state space to linearized space
    c=0.02567;
    x=[y(1);y(5);y(3);(y(2)+y(6)*sin(y(3)))/cos(y(3));y(6);y(4)/c];
end

