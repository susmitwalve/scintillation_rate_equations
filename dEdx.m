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
l=25;
imp=logspace(-2,0,l);
Einc1=imp;
dx=zeros(1,l);
Einc=imp;
dxf=0.5e-6;
for j=1:l
    if imp(j)<=1e-1
        dx(j)=1e-6;
    elseif imp(j)>1e-1 && imp(j)<=3e-1
        dx(j)=1e-4;
    else
        dx(j)=1e-3;
    end
end
for k = 1:length(Einc1)
    disp(k)
    E=array(:,1);
    dedx=array(:,2);
    i=0;
    while i<50000 && Einc(k)>0
        i=i+1;
        ddxcur=interp1(E,dedx,Einc(k),'spline');
        ddx(i)=ddxcur;
        x(i)=i*dx(k);
        Efin(i)=Einc(k);
        Einc(k)=Einc(k)-dx(k)*ddxcur;
    end
    ddx(i+1)=(Einc(k)+dx(k)*ddxcur)/dx(k);
    ddx(i+2)=0;
    x(i+1)=(i+1)*dx(k);
    x(i+2)=(i+2)*dx(k);
    m=1;
    xf(1)=1.01*dx(k);
    while xf(m)<x(length(x))
        m=m+1;
        xf(m)=xf(1)+(m-1)*dxf;
    end
    ddxf=interp1(x,ddx,xf(1:m-1));
    ddxf(m)=0;
    narr=(4.51*ddxf*1e6)/(pi*((3e-7)^2)*8.9);
    LY1(k)=0;
    for i=1:length(ddxf)
        ddxt=[1,2,3,5,7,9,11,15,18,22,24];
        nat=(4.51*ddxt*1e6)/(pi*((3e-7)^2)*8.9);
        LY=[1.333932236937050e-06,1.473029052131873e-06,1.554696141879172e-06,1.658606233702162e-06,1.724249037157726e-06,1.768142216255179e-06, 1.797968583870287e-06,1.831328677336973e-06, 1.841900262978692e-06,1.844270832667142e-06,1.841960838228944e-06];
        LY1(k)=LY1(k)+((dxf)*interp1(nat,LY,narr(i),'spline')/(length(ddxf)));
    end
    clear('ddx')
    clear('x')
    clear('Efin')
    clear('ddxf')
    clear('xf')
end
val=interp1(Einc1*1000,LY1,200);
semilogx(Einc1*1000,LY1/val,'*')
axis([1 1000 0.4 1.5])
grid on