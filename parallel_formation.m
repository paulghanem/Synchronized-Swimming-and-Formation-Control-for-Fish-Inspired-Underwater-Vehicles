clear all;
close all;
clc

addpath('./export_fig');

% parameters
d = 1;                % drag
m = 1;                  % mass
l = 1;                  % length
J = 1;                  % moment of inertia
Jd = J;                 % moment of inertia (disk)
%K = [.5 1 5 0];              % control gains
K = [.5 1 5];              % AW: control gain has extra zero?

% desired turning rate
b = J+m*l^2;            % paralell axis theorem
N=8;
% integrate

% AW : fix the seed for repatable results
rng('default');
rng(3);

thetad =45/180*pi*ones(N,1);
thetaddot=0*ones(N,1);

P=diag(N*ones(1,N))-1*(ones(1,N)'*ones(1,N));

omega=0*rand(N,1);
v=0*rand(N,1);
theta=pi*rand(N,1);
x=10*rand(N,1);
y=10*rand(N,1);
z0=[v,omega,theta,thetad,x,y];



ms = 12; % marker size
fs = 18; % font size


%tspan = [0 50];
%options = odeset('reltol',1e-3,'abstol',1e-6);

% AW : specify for points in tspan for smoother curves
tspan = linspace(0,50,1e4);
options = odeset('reltol',1e-3,'abstol',1e-6);

tic
[t,z] = ode45(@myodefun,tspan,z0,[],d,m,l,b,K,thetad,thetaddot,P,N);
toc


figure(1)
clf
plot(z(:,1:N),z(:,N+1:2*N))
hold on;
plot(z(end,1:N),z(end,N+1:2*N),'kx','MarkerSize',ms,'linewidth',2)
xlabel('Speed, v (m/s)')
ylabel('Angular Rate, \omega (rad/s)')
set(gca,'FontSize',fs)
set(gcf,'Color','w');
axis tight;
xlim([0 1.5]);
ylim([-1.5 1.5]);
axis square;
%ylim([2 2])
%export_fig v_omega.pdf


figure(2)
clf
plot(z(:,2*N+1:3*N),z(:,N+1:2*N))
set(gca,'FontSize',fs)
hold on
xlim([-1 2.5])
ylim([-1.5 1.5])
%xticks([-1:1:2.5]);
%yticks([-1.5:1:1.5]);

%axis equal;
plot(z(end,2*N+1:3*N),z(end,N+1:2*N),'kx','MarkerSize',ms,'linewidth',2)
%axis([-pi pi -2*pi 2*pi])
xlabel('Heading Angle, \theta (rad)')
ylabel('Angular Rate, \omega (rad/s)')
set(gca,'FontSize',fs)
%grid on
set(gcf,'Color','w');
axis square;
%ylim([2 2])
%export_fig theta_omega.pdf


figure(3)
clf
plot(t,z(:,2*N+1:3*N),t,z(:,3*N+1:4*N))
xlabel('Time (s)')
ylabel('Heading Angle, \theta (rad)')
set(gca,'FontSize',fs)
set(gcf,'Color','w');
%ylim([2 2])
xlim([tspan(1) tspan(end)])
axis square;
%export_fig t_theta.pdf


% figure(4)
% plot(t,abs(sum((exp(i*z(:,N+1:2*N)))/N,2)))
% %  figure(5)
% %  plot(t,mean(z(:,1:N).*z(:,N+1:2*N))*ones(1,length(t)))
% %


figure(5), clf,

plot(z(:,4*N+1:5*N),z(:,5*N+1:6*N))
xlabel('X Position (m)')
ylabel('Y Position (m)')
set(gca,'FontSize',fs)
xlim([0 30]);
ylim([0 30]);
xticks([0:5:30]);
yticks([0:5:30]);
axis square;
set(gcf,'Color','w');
%ylim([2 2])
export_fig x_y.pdf

% hold on
% plot((mean(z(:,N+1:2*N)))*cos(mean(z(:,N+1:2*N))*t),(mean(z(:,N+1:2*N)))*sin((mean(z(:,N+1:2*N))*t)))
% hold off

%%%%%%%%%%
% z = [v omega theta alphadot alpha x y thetad]
function zdot = myodefun(t,z,d,m,l,b,K,thetad,thetaddot,P,N)





% v dot 
zdot(1:N,1) = l*(z(N+1:2*N)).^2-d*z(1:N);
T(1:N,1) = (-K(1)*b*(z(N+1:2*N)))+b*K(2)*sin(z(2*N+1:3*N)-z(3*N+1:4*N))+(K(3)/N)*P*(z(N+1:2*N));

% omega dot
zdot(N+1:2*N,1) =-m*l/b*z(1:N).*(z(N+1:2*N))-T/b;

% theta dot
zdot(2*N+1:3*N,1) = z(N+1:2*N);

zdot(3*N+1:4*N,1)=thetaddot;

% dubins car
% zdot(4*N+1:5*N,1) = z(1:N,1).*cos(z(2*N+1:3*N));
% zdot(5*N+1:6*N,1) = z(1:N,1).*sin(z(2*N+1:3*N));

% Chaplygin sleigh (Jinseong paper)
omega = z(N+1:2*N);
theta = z(2*N+1:3*N);
zdot(4*N+1:5*N,1) = z(1:N,1).*cos(z(2*N+1:3*N)) - l*omega.*sin(theta);
zdot(5*N+1:6*N,1) = z(1:N,1).*sin(z(2*N+1:3*N)) + l*omega.*cos(theta);

end



