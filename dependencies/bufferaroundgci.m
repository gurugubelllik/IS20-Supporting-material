function X=bufferaroundgci(s,gci,nl,nr)

%%% Usage: X = bufferaroundgci(s,gci,nl,nr)

    gci=gci(:);
    Ns=length(s);
    win=[-nl:nr];
    win=win(:)';
    ndx=gci*ones(1, nl+nr+1) + win(ones(length(gci),1),:);
    ndx(ndx<1)=1;
    ndx(ndx>Ns)=Ns;
    X=s(ndx);

end