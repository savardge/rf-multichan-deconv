function [y]=bpfilt(x,dt,lf,hf)
%  BPFILT Bandpass filter time series. 
%         BPFILT(X,DT,LF,HF) takes a time series sampled at DT 
%	  and filters it with a 2nd order, 2 pass butterworth 
%         filter between frequencies LF and HF. If X is a matrix
%         BPFILT filters the individual rows of X.

% row vector
x=x(:)';

nyq=0.5/dt;
wn=[lf/nyq,hf/nyq];
[b,a]=butter(2,wn);
nx=size(x,1);
y=zeros(size(x));
for ix=1:nx
    y(ix,:)=filtfilt(b,a,x(ix,:));
end

end