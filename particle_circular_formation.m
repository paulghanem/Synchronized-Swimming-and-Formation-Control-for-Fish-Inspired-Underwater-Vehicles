clear all
clc
N=12;
M=N;
Tend=50;
Ts=1/5;
r=[1 0;0 2;-1 4; 3.5 0.5 ; -1 -5 ;1 -2; 0.5 0.5 ; 5 2; 4 8; 1 2; -3 -2; 0 0];
q=[0 pi/4 pi/2 2*pi/3 pi -pi/2 -pi/3 -pi/4 -2*pi/3 pi/6 -pi/6 5*pi/6]';
qd=zeros(N,1);

%K=0.2;
K0=2;


w0=0.3;


for tk=1:1:Tend/Ts
    R=sum(r,1)/N;
    
rb= r- R;

    for k=1:N
            for j=1:N
u(j,k)=(K0/N)*sin(q(j)-q(k));
            end 
    end 
    rd = [cos(q) , sin(q)];
    r=rd*Ts+r;
    X(:,tk)=r(:,1);
    Y(:,tk)=r(:,2);
    qd=sum(u,1);
    qd=qd';
   
    


q=qd*Ts+q;

end 

       for i= 1:N
   plot(X(i,:),Y(i,:));
   hold on 
    quiver(X(i,tk),Y(i,tk),rd(i,1)*5,rd(i,2)*5,'marker','o','color','k','Autoscale','off')
       end 
       xlabel('X (m)');
        ylabel('Y (m)');
        set(gca,'FontSize',20);
hold off
