%% Transformation Check
% This code checks whether the transformation obtained using 
% the dynamic extension is indeed diffeomorphic.

clear
clc
close all

% c=0.02567;
% a1=-0.28;a2=-0.993;a3=-0.1729;
% b1=0.0538;b2=-0.0122;b3=-0.0349;

syms x1 x2 x3 x4 x5 x6 a1 a2 a3 b1 b2 b3 c

z=[x1,
    x4*cos(x3)-x5*sin(x3),
    (a1*x4+(b1-c)*x5*x6)*cos(x3)-(a2*x5+(b2+c)*x4*x6)*sin(x3),
    x2,
    x4*sin(x3)+x5*cos(x3),
    (a1*x4+(b1-c)*x5*x6)*sin(x3)+(a2*x5+(b2+c)*x4*x6)*cos(x3)];

%jacobian(z,[x1,x2,x3,x4,x5,x6])
D=det(jacobian(z,[x1,x2,x3,x4,x5,x6]))