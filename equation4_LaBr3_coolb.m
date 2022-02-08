function res=equation4_LaBr3_coolb(ti,NSTE,NCe,t,constRCe,constS1,constK2Ce)
res=constS1*NSTE(ti,:)-constRCe*NCe(ti,:)-(constK2Ce*NCe(ti,:).*NCe(ti,:)/((t)^(0.5)));
end