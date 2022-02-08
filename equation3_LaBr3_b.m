function res=equation3_LaBr3_b(ti,ne,nh,NSTE,B,S1,K2E,RSTE,QSTE)
res=(B*ne(ti,:).*nh(ti,:))-((S1+RSTE+QSTE)*NSTE(ti,:))-(K2E*NSTE(ti,:).*NSTE(ti,:));
end