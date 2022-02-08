clear;
% RUN THE FILE reset_file.mat BEFORE EXECUTING THIS CODE
C=Constants2;
format short
load('index_file')
LY_mat=readmatrix('LY_file.txt');
ddx=[1,2,3,5,7,9,11,15,18,22,24];
n0_mat=(5.29*ddx*1e6)/(pi*((4*1e-7)^2)*13.7);
tx=4800;
num=0;
endx=600;
dx=(2*endx/(tx-1))*1e-7;
x2=linspace(-endx,endx,tx);
x222=abs(x2);
x222sp2=x222(1,(tx/2)+1:end).^2-x222(1,(tx/2):end-1).^2;
t=0;
ne=zeros(2,tx);
nh=zeros(2,tx);
NSTE=zeros(2,tx);
NCe=zeros(2,tx);
E_mat=zeros(tx/2,tx/2);
n0=n0_mat(index);
for d0=(tx/2)+1:tx
    ne(1,d0)=0.4*nborn(n0,x2(d0),C.rtrack);
    nh(1,d0)=0.4*nborn(n0,x2(d0),C.rtrack);
    NSTE(1,d0)=0.6*nborn(n0,x2(d0),C.rtrack);
end

for i=(tx/2)+1:tx
    for j=(tx/2)+1:tx
        if j==i-1
            E_mat(i-(tx/2),j-(tx/2))=-1/2;
        elseif j==i+1
            E_mat(i-(tx/2),j-(tx/2))=1/2;
        elseif j==i
            E_mat(i-(tx/2),j-(tx/2))=(dx/(x222(i)*1e-7));
        end
    end
end
E_mat(1,2)=1;
E_mat(1,1)=0;
E_inv=inv(E_mat);
ti=1;

Efield=zeros(1,tx);
e=1.60217662e-19;
dt=0.25e-16;
De=C.Dehot;
i1=0;
decayt=0;
mue=0;
Dh=C.Dh;
muh=0;
DSTE=C.DSTE;
B=0;
K3=C.K3;
RCe=C.RCe;
RSTE=C.RSTE;
QSTE=C.QSTE;
K2E=C.K2E;
K2Ce=C.K2Ce;
S1=0.72*C.Sfast+0.28*C.Sslow;
% S1=0;
j1=0;
while t<1e-12
    i1=i1+1;
    rho=(nh(ti,:)-ne(ti,:))*1e6*(e*(1e-4))/(5*8.85418782e-12);
    Efield((tx/2)+1:tx)=Ecalculator(tx,dx,E_inv,rho);
    for dn=(tx/2)+1:tx-1
        ne(2,dn)=equation1_LaBr3_a(Efield,ti,dn,dx,ne,tx,x222,De,mue);
        nh(2,dn)=equation2_LaBr3_a(Efield,ti,dn,dx,nh,tx,x222,Dh,muh);
        NSTE(2,dn)=equation3_LaBr3_a(ti,dn,dx,NSTE,x222,DSTE,tx);
    end
    ne(2,:)=ne(1,:)+(dt*(ne(2,:)+equation1_LaBr3_b(ti,ne,nh,B,K3)));
    nh(2,:)=nh(1,:)+(dt*(nh(2,:)+equation2_LaBr3_b(ti,ne,nh,B,K3)));
    NSTE(2,:)=NSTE(1,:)+(dt*(NSTE(2,:)+equation3_LaBr3_b(ti,ne,nh,NSTE,B,S1,K2E/((t+1e-16)^(0.5)),RSTE,QSTE)));
    NCe(2,:)=NCe(1,:)+(dt*(equation4_LaBr3_b(ti,NSTE,NCe,RCe,S1,K2Ce/((t+1e-16)^(0.5)))));
%     num=num+dt*R1E*2*pi*dx*sum(x222.*N(2,:));
    num=num+dt*RCe*2*pi*dx*sum(x222.*NCe(2,:));
    if mod(i1,50)==0
        j1=j1+1;
        decayt(1,j1)=RCe*dt*sum(x222sp2.*mean([NCe(2,(tx/2):end-1);NCe(2,(tx/2)+1:end)]))/(x222(tx)^2);
        decayt(2,j1)=t;
    end
    t=t+dt;    
    ne(1,:)=ne(2,:);
    nh(1,:)=nh(2,:);
    NSTE(1,:)=NSTE(2,:);
    NCe(1,:)=NCe(2,:);
    disp(t)
end
e=1.60217662e-19;

mue=C.mue;
De=C.De;
Dh=C.Dh;
muh=C.muh;
DSTE=C.DSTE;
B=C.B;
K3=C.K3;
RCe=C.RCe;
K2E=C.K2E;
K2Ce=C.K2Ce;
dtmax=0.98*((dx^2)/(2*De));
S1=0.72*C.Sfast+0.28*C.Sslow;
while t<0.1e-9
    i1=i1+1;
    rho=(nh(ti,:)-ne(ti,:))*1e6*(e*(1e-4))/(5*8.85418782e-12);
    Efield((tx/2)+1:tx)=Ecalculator(tx,dx,E_inv,rho);
    for dn=(tx/2)+1:tx-1
        ne(2,dn)=equation1_LaBr3_a(Efield,ti,dn,dx,ne,tx,x222,De,mue);
        nh(2,dn)=equation2_LaBr3_a(Efield,ti,dn,dx,nh,tx,x222,Dh,muh);
        NSTE(2,dn)=equation3_LaBr3_a(ti,dn,dx,NSTE,x222,DSTE,tx);
    end
    ne(2,:)=ne(1,:)+(dt*(ne(2,:)+equation1_LaBr3_b(ti,ne,nh,B,K3)));
    nh(2,:)=nh(1,:)+(dt*(nh(2,:)+equation2_LaBr3_b(ti,ne,nh,B,K3)));
    NSTE(2,:)=NSTE(1,:)+(dt*(NSTE(2,:)+equation3_LaBr3_b(ti,ne,nh,NSTE,B,S1,K2E/((t+1e-16)^(0.5)),RSTE,QSTE)));
    NCe(2,:)=NCe(1,:)+(dt*(equation4_LaBr3_b(ti,NSTE,NCe,RCe,S1,K2Ce/((t+1e-16)^(0.5)))));
    num=num+dt*RCe*2*pi*dx*sum(x222.*NCe(2,:));
    if mod(i1,50)==0
        j1=j1+1;
        decayt(1,j1)=RCe*dt*sum(x222sp2.*mean([NCe(2,(tx/2):end-1);NCe(2,(tx/2)+1:end)]))/(x222(tx)^2);
        decayt(2,j1)=t;
    end
    dt=dt_set(ne,nh,NSTE,NCe,dt,dtmax,tx);
    t=t+dt;    
    ne(1,:)=ne(2,:);
    nh(1,:)=nh(2,:);
    NSTE(1,:)=NSTE(2,:);
    NCe(1,:)=NCe(2,:);
    disp(t)
end
ne=zeros(size(ne));
dtmax=0.49*((dx^2)/(2*DSTE));

while t<1e-6
    i1=i1+1;
    rho=(nh(ti,:)-ne(ti,:))*1e6*(e*(1e-4))/(5*8.85418782e-12);
    Efield((tx/2)+1:tx)=Ecalculator(tx,dx,E_inv,rho);
    for dn=(tx/2)+1:tx-1
        nh(2,dn)=equation2_LaBr3_a(Efield,ti,dn,dx,nh,tx,x222,Dh,muh);
        NSTE(2,dn)=equation3_LaBr3_a(ti,dn,dx,NSTE,x222,DSTE,tx);
    end
    nh(2,:)=nh(1,:)+(dt*(nh(2,:)+equation2_LaBr3_b(ti,ne,nh,B,K3)));
    NSTE(2,:)=NSTE(1,:)+(dt*(NSTE(2,:)+equation3_LaBr3_b(ti,ne,nh,NSTE,B,S1,K2E/((t+1e-16)^(0.5)),RSTE,QSTE)));
    NCe(2,:)=NCe(1,:)+(dt*(equation4_LaBr3_b(ti,NSTE,NCe,RCe,S1,K2Ce/((t+1e-16)^(0.5)))));  
    num=num+dt*RCe*2*pi*dx*sum(x222.*NCe(2,:));
    if mod(i1,5)==0
        j1=j1+1;
        decayt(1,j1)=RCe*dt*sum(x222sp2.*mean([NCe(2,(tx/2):end-1);NCe(2,(tx/2)+1:end)]))/(x222(tx)^2);
        decayt(2,j1)=t;
    end
    dt=dt_set_n0(nh,NSTE,NCe,dt,dtmax,tx);
    t=t+dt;    
    nh(1,:)=nh(2,:);
    NSTE(1,:)=NSTE(2,:);
    NCe(1,:)=NCe(2,:);
    disp(t)
end

LY_mat(index)=num;
disp("LY: "+num2str(num))
if index==1
    writematrix(decayt,'decayt1.txt')
elseif index==2
    writematrix(decayt,'decayt2.txt')
elseif index==3
    writematrix(decayt,'decayt3.txt')
elseif index==4
    writematrix(decayt,'decayt5.txt')
elseif index==5
    writematrix(decayt,'decayt7.txt')
elseif index==6
    writematrix(decayt,'decayt9.txt')
elseif index==7
    writematrix(decayt,'decayt11.txt')
elseif index==8
    writematrix(decayt,'decayt15.txt')
elseif index==9
    writematrix(decayt,'decayt18.txt')
elseif index==10
    writematrix(decayt,'decayt22.txt')
elseif index==11
    writematrix(decayt,'decayt24.txt')
end
index=index+1;
save('index_file','index')
writematrix(LY_mat,'LY_file.txt')

scintillation_LaBr3

function y=nborn(ninc,x,r0)
y=ninc*exp(-((x/r0)^2));
end
