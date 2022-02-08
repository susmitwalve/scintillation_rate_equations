function res=equation1_LaBr3_b(ti,ne,nh,B,K3)
res=-B*ne(ti,:).*nh(ti,:)-K3*ne(ti,:).*ne(ti,:).*nh(ti,:);
end