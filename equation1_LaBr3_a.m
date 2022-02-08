function res1=equation1_LaBr3_a(Efield,ti,xi,dx,ne,tx,x222,De,mue)
res1=0;
if xi==(tx/2)+1
    res1=res1+De*(((ne(ti,xi+2)+ne(ti,xi)-2*ne(ti,xi+1))/dx^2)+((ne(ti,xi+1)-ne(ti,xi))/(x222(xi)*1e-7*dx)))+mue*((ne(ti,xi+1)*Efield(xi+1)*x222(xi+1))-(ne(ti,xi)*Efield(xi)*x222(xi)))/(x222(xi)*dx);
else
    res1=res1+De*((((2*x222(xi)*1e-7)+dx)*ne(ti,xi+1))+(((2*x222(xi)*1e-7)-dx)*ne(ti,xi-1))-(4*x222(xi)*1e-7*ne(ti,xi)))/(2*x222(xi)*1e-7*(dx^2))+mue*((ne(ti,xi+1)*Efield(xi+1)*x222(xi+1))-(ne(ti,xi-1)*Efield(xi-1)*x222(xi-1)))/(2*x222(xi)*dx);
end
end