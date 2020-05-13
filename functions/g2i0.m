function [ip,ib]=g2i0(nib,np)

% Function takes grid index and produces corresponding
% slowness and backazimuthal indices.
ib=fix((nib-1)/np)+1;
ip=nib-(ib-1)*np;

return

