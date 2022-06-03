function [circleHandle,z0] = circle(x0,y0,r0,n)
% [circleHandle,z0] = circle(x0,y0,r0,n) function plots a circle of radius r0 at center (x0,y0). circle perimeter
% has 100 points unless specified by number n
circleHandle = nan;
if nargin<4
    n = 100;
end
theta = linspace(0,2*pi,n);
z =r0*exp(1j*theta); % z holds complex x and y points of circle at origin
z0 = z +x0+1j*y0; % shift the circle to where it's supposed to be
if nargout<2
    circleHandle = plot(z0);
end
end