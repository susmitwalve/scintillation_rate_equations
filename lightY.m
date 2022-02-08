function LY1=lightY(dxf,len,dect,ddxt,decaymat,ddxf,LY1)
for k=1:length(ddxf)%at a particular value of dE/dx
    for jt=1:len%at a particular time
        dect(jt)=interp1(ddxt,decaymat(:,jt),ddxf(k),'spline');
        
    end
    LY1=LY1+((dxf)*dect/(length(ddxf)));
    disp(k)
end
end