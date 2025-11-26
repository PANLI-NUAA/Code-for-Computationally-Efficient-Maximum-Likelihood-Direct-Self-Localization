function [a,b]=cosine_surrogate(c,theta0)
a=abs(c)/2;
b=-abs(c)*theta0-c*sin(theta0);