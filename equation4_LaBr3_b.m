function res=equation4_LaBr3_b(ti,NSTE,NCe,RCe,S1,K2Ce)
res=S1*NSTE(ti,:)-RCe*NCe(ti,:)-(K2Ce*NCe(ti,:).*NCe(ti,:));
end