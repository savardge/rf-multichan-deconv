function [y]=taper(x,nt,dt,t1,t2)
%   TAPER  Taper time series.
%          TAPER(X,NT,DT,T1,T2) takes time series sampled at DT
%          and tapers it with a cosine taper NT seconds long from
%          begin point T1-NT and with reverse cosine taper from
%          point T2 to point T2+NT. Points outside the range
%          (T1-NT,T2+NT) are zeroed. If T1/T2 is negative then
%          taper is not implemented at beginning/end. If X is an
%          array of seismograms, then the taper is applied to each
%          row of X.
%
nx=max(size(x));
ns=min(size(x));
taper=ones(1,nx);
it=[0:fix(nt/dt)]*dt/nt;
ct=0.5-0.5*cos(pi*it);
it1=fix(t1/dt+1);
it2=fix(t2/dt+1);
if t1 > 0
    taper(it1-fix(nt/dt):it1)=ct;
    taper(1:it1-fix(nt/dt))=zeros(size(1:it1-fix(nt/dt)));
end
if t2 > 0
    if t2 > nx*dt-nt
        [t2, nx*dt-nt]
        t2 - (nx*dt-nt)
        error ('T2 > NX * DT - NT')
    end
    taper(it2:it2+fix(nt/dt))=fliplr(ct);
    taper(it2+fix(nt/dt):nx)=zeros(size(taper(it2+fix(nt/dt):nx)));
end
for is=1:ns
    y(is,:)=taper.*x(is,:);
end
