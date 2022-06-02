clear all;
close all;
clc

rng(9)
% parameters
d = 1;                % drag
m = 1;                  % mass
l = 1;                  % length
J = 1;                  % moment of inertia
Jd = J;                 % moment of inertia (disk)
K = [-.5 4 0 0];              % control gains
% desired turning rate
b = J+m*l^2;            % paralell axis theorem
N=8;
% integrate

thetad =45/180*pi*ones(N,1);
thetaddot=0*ones(N,1);

P=diag(N*ones(1,N))-1*(ones(1,N)'*ones(1,N));

omega=0*rand(N,1);
v=0*rand(N,1);
theta=2*pi*rand(N,1);
x=10*rand(N,1);
y=10*rand(N,1);
%x(7) = x(7)-2;
% move 7th fish to prevent overlap
z0=[v,omega,theta,thetad,x,y];




%tspan = [0 50];
tspan = linspace(0,35,1e4);
options = odeset('reltol',1e-3,'abstol',1e-6);
tic
[t,z] = ode45(@myodefun,tspan,z0,[],d,m,l,b,K,thetad,thetaddot,P,N);
toc


figure(1)
clf
plot(z(:,1:N),z(:,N+1:2*N))
xlabel('v (m/s)')
ylabel('\omega (rad/s)')


figure(2)
clf
plot(z(:,2*N+1:3*N),z(:,N+1:2*N))

hold on


plot(z(end,2*N+1:3*N),z(end,N+1:2*N),'kx')

axis([-pi pi -2*pi 2*pi])
xlabel('\theta (rad)')
ylabel('\omega (rad/s)')
grid on
axis square

figure(3)
clf

plot(t,z(:,2*N+1:3*N),t,z(:,3*N+1:4*N))
xlabel('t (s)')
ylabel('\theta (rad)')

figure(4)
plot(t,abs(sum((exp(i*z(:,N+1:2*N)))/N,2)))
%  figure(5)
%  plot(t,mean(z(:,1:N).*z(:,N+1:2*N))*ones(1,length(t)))
%
figure(5), clf,

plot(z(:,4*N+1:5*N),z(:,5*N+1:6*N))
xlabel('position x (m)')
ylabel('position y (m)')
% hold on
% plot((mean(z(:,N+1:2*N)))*cos(mean(z(:,N+1:2*N))*t),(mean(z(:,N+1:2*N)))*sin((mean(z(:,N+1:2*N))*t)))
% hold off



% rename variables
v = z(:,1:N);
w = z(:,N+1:2*N);
theta = z(:,2*N+1:3*N);
theta_d = z(:,3*N+1:4*N);
x = z(:,4*N+1:5*N);
y = z(:,5*N+1:6*N);

% get fish coordinates (from Brian Free)
% draw fish
tj = .02; % twice the maximum thickness from camberline, m
hj = .0;
L = 0.15;
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
xFish = fliplr(real(zFish));
yFish = fliplr(imag(zFish));

% movie
movieFlag = 1;
numPts = length(x);
hold on 



skip = 40;
c = lines(8);    
    for i = numPts
        % plot planar motion
       
        colors = get(gca,'colororder');
        % current position
        plot( x(i,:) , y(i,:), 'k.', 'linewidth',2,'MarkerSize',20); hold on;
       
        % trail
        trailSize = 1;
        nStart = max(1, i-trailSize);
        if ( i > 1 )
        plot( x(nStart:i,:) , y(nStart:i,:), 'linewidth',3);        
        end
        % current pose
        L = 0.31;
        for j = 1:1:N
            xc = [x(i,j) y(i,j)];
            th = pi + theta(i,j); % add pi since fish is facing negative x axis
            R = [cos(th) -sin(th); sin(th) cos(th)];
            scale = 10;
            for k = 1:1:length(xFish)
                shapeFish = scale*R*[xFish(k) 3*yFish(k)]';                
                xbody(k) = xc(1) + shapeFish(1);
                ybody(k) = xc(2) + shapeFish(2);
            end
            %xbody = [xc(1)-L*cos(th) xc(1)+L*cos(th)];
            %ybody = [xc(2)-L*sin(th) xc(2)+L*sin(th)];
            fill( xbody, ybody, c(j,:), 'linewidth',1); hold on;
        end         
        xlabel('X (m)');
        ylabel('Y (m)');
        set(gca,'FontSize',16);
        axis equal;
        xlim([-15 10]);
        ylim([-15 10]);          
        hold off;
        
        % 
        
        
        drawnow;
       
             
    end
    




%%%%%%%%%%
% z = [v omega theta alphadot alpha x y thetad]
function zdot = myodefun(t,z,d,m,l,b,K,thetad,thetaddot,P,N)






r=[z(4*N+1:5*N,1),z(5*N+1:6*N,1)];

  R= 1/N*[ sum(z(4*N+1:5*N,1)), sum(z(5*N+1:6*N,1)) ];    
   rb=r-R ;       

for k=1:N
            for j=1:N
u(j,k)=sin(z(2*N+j)-z(2*N+k));
            end 
end 
    qd=sum(u,1);
rd=[z(1:N).*cos(z(2*N+1:3*N)),z(1:N).*sin(z(2*N+1:3*N))];
for j=1:N
           rs(j)=rb(j,:)*rd(j,:)';
     end
zdot(1:N,1) = l*(z(N+1:2*N)).^2-d*z(1:N);
T(1:N,1) = K(1)*b*(z(N+1:2*N)) -K(2)*b*qd(:)+(K(3)/N)*P*(z(N+1:2*N));
  
zdot(N+1:2*N,1) =-m*l/b*z(1:N).*((z(N+1:2*N)))-T/b;
zdot(2*N+1:3*N,1) = z(N+1:2*N);



zdot(3*N+1:4*N,1)=thetaddot;

zdot(4*N+1:5*N,1) = z(1:N,1).*cos(z(2*N+1:3*N));
zdot(5*N+1:6*N,1) = z(1:N,1).*sin(z(2*N+1:3*N));
end





