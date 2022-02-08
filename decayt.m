clear;
array=[1.000E-03, 3.167E+01, 
1.500E-03, 3.025E+01, 
2.000E-03, 2.756E+01, 
2.500E-03, 2.510E+01, 
3.000E-03, 2.300E+01, 
3.500E-03, 2.123E+01, 
4.000E-03, 1.973E+01, 
4.500E-03, 1.845E+01, 
5.000E-03, 1.734E+01, 
6.000E-03, 1.553E+01, 
7.000E-03, 1.410E+01, 
8.000E-03, 1.294E+01, 
9.000E-03, 1.199E+01, 
1.000E-02, 1.118E+01, 
1.500E-02, 8.499E+00, 
2.000E-02, 6.970E+00, 
2.500E-02, 5.972E+00, 
3.000E-02, 5.265E+00, 
3.500E-02, 4.736E+00, 
4.000E-02, 4.324E+00, 
4.500E-02, 3.993E+00, 
5.000E-02, 3.721E+00, 
6.000E-02, 3.302E+00, 
7.000E-02, 2.991E+00, 
8.000E-02, 2.753E+00, 
9.000E-02, 2.563E+00, 
1.000E-01, 2.409E+00, 
1.500E-01, 1.934E+00, 
2.000E-01, 1.692E+00, 
2.500E-01, 1.548E+00, 
3.000E-01, 1.454E+00, 
3.500E-01, 1.390E+00, 
4.000E-01, 1.345E+00, 
4.500E-01, 1.312E+00, 
5.000E-01, 1.287E+00, 
6.000E-01, 1.256E+00, 
7.000E-01, 1.240E+00, 
8.000E-01, 1.232E+00, 
9.000E-01, 1.230E+00, 
1.000E+00, 1.231E+00, 
1.500E+00, 1.266E+00, 
2.000E+00, 1.317E+00, 
2.500E+00, 1.370E+00, 
3.000E+00, 1.424E+00, 
3.500E+00, 1.476E+00, 
4.000E+00, 1.528E+00, 
4.500E+00, 1.579E+00, 
5.000E+00, 1.629E+00, 
6.000E+00, 1.728E+00, 
7.000E+00, 1.825E+00, 
8.000E+00, 1.921E+00, 
9.000E+00, 2.017E+00, 
1.000E+01, 2.112E+00, 
1.500E+01, 2.584E+00, 
2.000E+01, 3.056E+00, 
2.500E+01, 3.529E+00, 
3.000E+01, 4.005E+00, 
3.500E+01, 4.482E+00, 
4.000E+01, 4.961E+00, 
4.500E+01, 5.441E+00, 
5.000E+01, 5.924E+00, 
6.000E+01, 6.892E+00, 
7.000E+01, 7.864E+00, 
8.000E+01, 8.839E+00, 
9.000E+01, 9.817E+00, 
1.000E+02, 1.080E+01, 
1.500E+02, 1.572E+01, 
2.000E+02, 2.068E+01, 
2.500E+02, 2.565E+01, 
3.000E+02, 3.062E+01, 
3.500E+02, 3.561E+01, 
4.000E+02, 4.060E+01, 
4.500E+02, 4.560E+01, 
5.000E+02, 5.060E+01, 
6.000E+02, 6.061E+01, 
7.000E+02, 7.064E+01, 
8.000E+02, 8.066E+01, 
9.000E+02, 9.070E+01, 
1.000E+03, 1.007E+02, 
1.500E+03, 1.510E+02, 
2.000E+03, 2.013E+02, 
2.500E+03, 2.516E+02, 
3.000E+03, 3.019E+02, 
3.500E+03, 3.522E+02, 
4.000E+03, 4.026E+02, 
4.500E+03, 4.529E+02, 
5.000E+03, 5.032E+02, 
6.000E+03, 6.039E+02, 
7.000E+03, 7.047E+02, 
8.000E+03, 8.054E+02, 
9.000E+03, 9.061E+02, 
1.000E+04, 1.007E+03];
k=1;
l=25;
imp=logspace(-2,0,l);
Einc1=imp;
dx=2e-2;
Einc=2.5;
dxf=2e-2;
E=array(:,1);
dedx=array(:,2);
i=0;
while i<50000 && Einc>0
    i=i+1;
    ddxcur=interp1(E,dedx,Einc,'spline');
    ddx(i)=ddxcur;
    x(i)=i*dx;
    Efin(i)=Einc;
    Einc=Einc-dx*ddxcur;
end
ddx(i+1)=(Einc+dx*ddxcur)/dx;
ddx(i+2)=0;
x(i+1)=(i+1)*dx;
x(i+2)=(i+2)*dx;
m=1;
xf(1)=1.01*dx;
while xf(m)<x(length(x))
    m=m+1;
    xf(m)=xf(1)+(m-1)*dxf;
end
ddxf=interp1(x,ddx,xf(1:m-1));
ddxf(m)=0;
narr=(5.29*ddxf*1e6)/(pi*((4e-7)^2)*13.7);
LY1(k)=0;
decayt1=readmatrix('decayt1.txt');
decayt2=readmatrix('decayt2.txt');
decayt3=readmatrix('decayt3.txt');
decayt5=readmatrix('decayt5.txt');
decayt7=readmatrix('decayt7.txt');
decayt9=readmatrix('decayt9.txt');
decayt11=readmatrix('decayt11.txt');
decayt15=readmatrix('decayt15.txt');
decayt18=readmatrix('decayt18.txt');
decayt22=readmatrix('decayt22.txt');
decayt24=readmatrix('decayt24.txt');
t=decayt1(2,:);
len=length(t);
decaymat=zeros(11,len);
decaymat(1,:)=interp1(decayt1(2,:),decayt1(1,:),t,'spline');
decaymat(2,:)=interp1(decayt2(2,:),decayt2(1,:),t,'spline');
decaymat(3,:)=interp1(decayt3(2,:),decayt3(1,:),t,'spline');
decaymat(4,:)=interp1(decayt5(2,:),decayt5(1,:),t,'spline');
decaymat(5,:)=interp1(decayt7(2,:),decayt7(1,:),t,'spline');
decaymat(6,:)=interp1(decayt9(2,:),decayt9(1,:),t,'spline');
decaymat(7,:)=interp1(decayt11(2,:),decayt11(1,:),t,'spline');
decaymat(8,:)=interp1(decayt15(2,:),decayt15(1,:),t,'spline');
decaymat(9,:)=interp1(decayt18(2,:),decayt18(1,:),t,'spline');
decaymat(10,:)=interp1(decayt22(2,:),decayt2(1,:),t,'spline');
decaymat(11,:)=interp1(decayt24(2,:),decayt24(1,:),t,'spline');
ddxt=[1,2,3,5,7,9,11,15,18,22,24];
dect=zeros(1,len);
LY1=zeros(1,len);
LY1=lightY(dxf,len,dect,ddxt,decaymat,ddxf,LY1);