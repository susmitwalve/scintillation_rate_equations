function res=equation3_LaBr3_a(ti,xi,dx,NSTE,x222,DE,tx)
if xi==(tx/2)+1
    res=DE*(((NSTE(ti,xi+2)+NSTE(ti,xi)-2*NSTE(ti,xi+1))/dx^2)+((NSTE(ti,xi+1)-NSTE(ti,xi))/(x222(xi)*1e-7*dx)));
else
    res=DE*((((2*x222(xi)*1e-7)+dx)*NSTE(ti,xi+1))+(((2*x222(xi)*1e-7)-dx)*NSTE(ti,xi-1))-(4*x222(xi)*1e-7*NSTE(ti,xi)))/(2*x222(xi)*1e-7*(dx^2));
end
end