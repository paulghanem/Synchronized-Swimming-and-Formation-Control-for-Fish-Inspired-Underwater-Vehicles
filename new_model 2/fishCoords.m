function [xFish,yFish] = fishCoords(x, y, L, angle)
% get fish coordinates (from Brian Free)
% draw fish
tj = L*.008/0.06; % twice the maximum thickness from camberline, m
hj = .0;
%L = 0.06;
joukowski.T = tj/L; % thickness ratio
joukowski.H = hj/L; % camber ratio
joukowski.c = L/4;
joukowski.zeta0 = joukowski.c*(-4/(3*sqrt(3))*joukowski.T+1j*2*joukowski.H); % center of lifting cylinder in zeta plane
joukowski.r0 = L*(1/4+joukowski.T/(3*sqrt(3))); % radius of lifting cylinder in zeta plane, m
[~,zetaCircle] = circle(real(joukowski.zeta0),imag(joukowski.zeta0),joukowski.r0);  % close the figure created by circle command
zFish = zetaCircle+joukowski.c^2./(zetaCircle);
zfog = zFish;
zp(1) = -146.48/2+9.77+10.47j; % values from CAD of 3D printed fish
zp(2) = -146.48/2+58.62+12.21j;
zp(3) = -146.48/2+9.77-10.47j;
zp(4) = -146.48/2+58.62-12.21j;
s = 0.001*zp; % convert from millimeters to meters
sog = s;
xFish = fliplr(real(zFish)) - L/2;
yFish = fliplr(imag(zFish));

xc = [x y];
th = pi + angle; % add pi since fish is facing negative x axis
R = [cos(th) -sin(th); sin(th) cos(th)];

scale = 1;
for k = 1:1:length(xFish)
    shapeFish = scale*R*[xFish(k) 3*yFish(k)]';
    xbody(k) = xc(1) + shapeFish(1);
    ybody(k) = xc(2) + shapeFish(2);
end

xFish = xbody;
yFish = ybody;
%fill( xbody, ybody, 'b', 'linewidth',1); hold on;    

end