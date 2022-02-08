function res=equation3_LaBr3_coolb(ti,ne,nh,NSTE,t,constB,constS1,constK2E)
res=(constB*ne(ti,:).*nh(ti,:))-(constS1*NSTE(ti,:))-(constK2E*NSTE(ti,:).*NSTE(ti,:)/((t)^(0.5)));
end