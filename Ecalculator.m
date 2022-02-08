function resV=Ecalculator(tx,dx,E_inv,rho)
res_mat=rho((tx/2)+1:tx)*dx;
% res_mat=rho*dx;
resV=transpose(E_inv*transpose(res_mat));
end