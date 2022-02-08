function res1=equation2_LaBr3_a(Efield,ti,xi,dx,nh,tx,x222,Dh,muh)
res1=0;
if xi==(tx/2)+1
    res1=res1+Dh*(((nh(ti,xi+2)+nh(ti,xi)-2*nh(ti,xi+1))/dx^2)+((nh(ti,xi+1)-nh(ti,xi))/(x222(xi)*1e-7*dx)))-muh*((nh(ti,xi+1)*Efield(xi+1)*x222(xi+1))-(nh(ti,xi)*Efield(xi)*x222(xi)))/(x222(xi)*dx);
else
    res1=res1+Dh*((((2*x222(xi)*1e-7)+dx)*nh(ti,xi+1))+(((2*x222(xi)*1e-7)-dx)*nh(ti,xi-1))-(4*x222(xi)*1e-7*nh(ti,xi)))/(2*x222(xi)*1e-7*(dx^2))-muh*((nh(ti,xi+1)*Efield(xi+1)*x222(xi+1))-(nh(ti,xi-1)*Efield(xi-1)*x222(xi-1)))/(2*x222(xi)*dx);
end
end