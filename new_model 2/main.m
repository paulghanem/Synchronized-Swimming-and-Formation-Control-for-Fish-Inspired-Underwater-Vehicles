clear; close all;  clc;

addpath('./brewer');
addpath('./cmocean_v2.0/cmocean');
addpath('./export_fig/');

lastk = 200; % this many points are used to plot final limit cycle


% parameters
m = 1.4;
l = 0.31;
d = 0.5;
b = 0.1395;
a = (m*l^2)/(b*d);
N = 8;

T = 20*120; % total simulation time

% parallel gains
k1 = 0.5;
k2 = 2;

% circular gains
k1c = k1;
k2c = k2;
k3c = -.3;

vref = 0*b*k1/(m*l);
omegaref = 0.2;


% create communication topology

%% graph used by Laplacian controllers 
% %all-to-all
% D = eye(N,N)*(N-1);
% A = ones(N,N);
% for i =1:1:N
%     A(i,i) = 0;
% end
% L = D-A;
% rng(3);

% circulant
L = zeros(N,N);
for i = 1:1:N
   L(i,i) = 2;
   % aft neighbor
   if ( i == 1)
       L(i,end) = -1;
   else
       L(i,i-1) = -1;
   end
   % fore neighbor
   if ( i == N)
       L(i,1) = -1;
   else
      L(i,i+1) = -1; 
   end       
end

%% ploting
% plot colors
shift = 2;
%cmap = parula(N); 
%cmap = cmocean('ice',N);
cmap = cmocean('deep',N+shift);

% shift
cmap = cmap(shift:end,:);

%% initial conditions

% %Random heading and random angular rate
% theta0 = rand([N 1])*2*pi;
% omega0 = rand([N 1])*1;

% %Random heading and zero angular rate
% r0 = rand(N,1)*10 + rand(N,1)*10*1i;
% v0 = 0*rand(N,1);
% 
% % theta0 = [0 45 90]'*pi/180;
% % omega0 = 0.1*[-1 0.5 1]';
% theta0 = rand([N 1])*pi;
% omega0 = rand(N,1)*0;

% Equal heading and equal angular rate
% theta0 = [ones(N,1)*0];
% omega0 = [ones(N,1)*1]

% Equal heading and random angular rate
% theta0 = [ones(N,1)*0];
% omega0 = rand([N 1])*1;

% % Np fish anti-parallel and random angular rate
% Np  = 1;
% theta0 = [ones([Np 1])*0; ones([N-Np 1])*pi ];
% omega0 = rand([N 1])*1;

% % circular initial formation
% th_circ = [0:2*pi/(N):2*pi];
% th_circ = th_circ(1:end-1)';
% R = 5;
% r0 = R*cos(th_circ) + R*sin(th_circ)*1i;
% v0 = 0*rand(N,1);
% theta0 = rand([N 1])*2*pi;
% omega0 = 0*rand(N,1);

% line initial formation
Linit = 5;
dely = 0;
r0 = linspace(-Linit/2,Linit/2,N)' - 1i*dely;
v0 = 0*rand(N,1);
theta0 = rand([N 1])*2*pi;
omega0 = 0*rand(N,1);


% % from paul's code
% zfix =[         0         0    1.3110         0    8.7639    9.5789;
%     0         0    1.7552         0    8.9461    5.3317;
%     0         0    0.4410         0    0.8504    6.9188;
%     0         0    0.6224         0    0.3905    3.1552;
%     0         0    2.5156         0    1.6983    6.8650;
%     0         0    3.0419         0    8.7814    8.3463;
%     0         0    0.9847         0    0.9835    0.1829;
%     0         0    2.1750         0    4.2111    7.5014];
% x0 = zfix(:,5);
% y0 = zfix(:,6);
% r0 = x0 + 1i * y0;
% theta0 = zfix(:,3);
% v0 = zeros(N,1);
% omega0 = zeros(N,1);

z0 = [r0; theta0; v0; omega0];

%% simulation 
tspan = [0:0.1:T];
%options = odeset('reltol',1e-12,'abstol',1e-12);
options = odeset('reltol',1e-3,'abstol',1e-3);
[t,zhist] = ode45(@parallelSys,tspan,z0,[],m,l,b,d,k1,k2,N,L,a,vref,omegaref,k1c,k2c,k3c);


%% plots
x_hist = real(zhist(:,1:N));
y_hist = imag(zhist(:,1:N));
th_hist = zhist(:,N+1:2*N);
v_hist = zhist(:,2*N+1:3*N);
omega_hist = zhist(:,3*N+1:4*N);


grayRGB = [1 1 1]*0.5;
figure;

hold on;
Lfish = 2;


for i = 1:1:N
    plot(x_hist(:,i),y_hist(:,i),'Color',cmap(i,:),'linewidth',1);
         
end
for i = 1:1:N
   [xFish,yFish] = fishCoords(x_hist(end,i), y_hist(end,i), Lfish, th_hist(end,i));
   fill(xFish,yFish,cmap(i,:), 'linewidth',1); hold on;   
   plot(xFish,yFish, 'k', 'linewidth',1); hold on; 
   plot(x_hist(end-lastk:end,i),y_hist(end-lastk:end,i),'m','Linewidth',2);
end
plot(x_hist(1,:),y_hist(1,:),'ko','MarkerFaceColor','k','MarkerSize',3)
xlabel('X (m)');
ylabel('Y (m)');
box on;
axis square;
axis equal;
axis tight;
set(gcf,'Color','w')
set(gca,'FontSize',16,'FontName','Arial')
%export_fig('parallel_fish_tracks.pdf')

% figure;
% for i = 1:1:N
% plot(t,mod(th_hist(:,i),2*pi)*180/pi,'linewidth',2,'Color',cmap(i,:),'linewidth',2);
% hold on;
% end
% ylabel('Heading (deg.)');
% xlabel('Time (sec.)');
% set(gca,'FontSize',16,'FontName','Arial')
% ylim([0 270]);
% xlim([0 90]);
% yticks([0:90:360]);
% set(gcf,'Color','w')
% %export_fig('parallel_t_vs_theta.pdf')

figure;
for i = 1:1:N
plot(v_hist(:,i), omega_hist(:,i)*180/pi,'Color',cmap(i,:));
hold on;
plot(v_hist(end-lastk:end,i),omega_hist(end-lastk:end,i)*180/pi,'m','Linewidth',2); 
plot(v_hist(1,i),omega_hist(1,i)*180/pi,'ko','MarkerFaceColor','k','MarkerSize',3)

hold on;
end

set(gca,'FontSize',16,'FontName','Arial')
set(gcf,'Color','w');
xlabel('Speed (m/s)');
axis square;
ylabel('Angular Rate (deg/s)');
%export_fig('parallel_v_vs_omega.pdf')

% figure;
% for i = 1:1:N
% plot(mod(th_hist(:,i),2*pi)*180/pi, omega_hist(:,i)*180/pi,'Color',cmap(i,:));
% hold on;
% plot(mod(th_hist(end-lastk:end,1),2*pi)*180/pi,omega_hist(end-lastk:end,1)*180/pi,'k','Linewidth',2); 
% plot(mod(th_hist(1,i),2*pi)*180/pi,omega_hist(1,i)*180/pi,'ko','MarkerFaceColor','k','MarkerSize',3)
% end
% set(gca,'FontSize',16,'FontName','Arial')
% set(gcf,'Color','w');
% xlabel('Heading (deg)');
% xlim([0 360]);
% axis square;
% ylabel('Angular Rate (deg/s)');
% xlim([0 360]);
% %export_fig('parallel_th_vs_omega.pdf')


function xdot = parallelSys(t,x,m,l,b,d,k1,k2,N,L,a,vref,omegaref,k1c,k2c,k3c)

% unpack state
r = x(1:N);
theta = x(N+1:2*N);
v = x(2*N+1:3*N);
omega = x(3*N+1:4*N);

%% Parallel: Laplacian control
u = zeros(N,1);
couplingTerm1 = zeros(N,1);
for k = 1:1:N
    Lk = L(k,:);
    val1 = 1i * exp(1i * theta(k) );
    val2 = Lk * exp(1i * theta );    
    couplingTerm1(k) = complexIP(val1,val2);    
end
u_p_lap = -b*k1*omega + b*k2/N*couplingTerm1;

%% Parallel: All-to-all control
couplingTerm2 = zeros(N,1);
for k = 1:1:N
    for j = 1:1:N
        couplingTerm2(k) = couplingTerm2(k) - sin( theta(j) - theta(k) );
    end
end
u_p_all = -b*k1*omega + b*k2/N*couplingTerm2;

%% Circular: Laplacian
couplingTerm_c_lap = zeros(N,1);
coefficient = ( (l/d)*omega  - vref/omegaref);
c = r + vref/omegaref*1i*exp(1i*theta);
for k = 1:1:N
    Lk = L(k,:);    
    val1 = exp(1i * theta(k) );
    val2 = Lk * c; %exp(1i * c);    
    couplingTerm_c_lap(k) = complexIP(val1,val2);    
end

 u_c_lap = -b*k1c*omega - b*k2c*sin(theta - omegaref*t) ...
     - b*k3c*coefficient.*couplingTerm_c_lap;
 
% figure(10);
% plot(real(r),imag(r),'bo'); hold on;
% plot(real(c),imag(c),'r+');
% hold on;

%% switch control from parallel/circular etc. here
u = u_c_lap; % circular, laplacian
%u = u_p_lap; % parallel laplacian
%u = u_p_all; % prallel all-to-all


%% dynamics
%omegadot = -a*omega.^3 + -u/b; % simplified model
omegadot = -m*l/b*v.*omega + -u/b; % actual model

% state rate
xdot(1:N,1) = v.*exp(1i*theta);% r-dot trailing edge
%xdot(1:N,1) = v.*exp(1i*theta) + l*omega.*1i.*exp(1i*theta);% r-dot for CM
xdot(N+1:2*N,1) = omega;%th-dot
xdot(2*N+1:3*N,1) = l*omega.^2 - d*v;%v-dot
xdot(3*N+1:4*N,1) = omegadot;

end



