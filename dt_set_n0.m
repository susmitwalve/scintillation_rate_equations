function dtnew=dt_set_n0(nh,NSTE,NCe,dtold,dtmax,tx)
nhdiff=max(abs(((nh(2,[(tx/2)+1,tx])/dtold-nh(2,[(tx/2)+1,tx])/dtold))./nh(2,[(tx/2)+1,tx])));
NSTEdiff=max(abs(((NSTE(2,[(tx/2)+1,tx])/dtold-NSTE(2,[(tx/2)+1,tx])/dtold))./NSTE(2,[(tx/2)+1,tx])));
NCediff=max(abs(((NCe(2,[(tx/2)+1,tx])/dtold-NCe(2,[(tx/2)+1,tx])/dtold))./NCe(2,[(tx/2)+1,tx])));
maxall=max([nhdiff,NSTEdiff,NCediff]);

dt1=0.1/maxall;

if dt1>dtmax
    dtnew=dtmax;
elseif isnan(dt1)
    dtnew=dtold;
else
    dtnew=dt1;
end
end