% Simple rectangular grid in ray parameter
% and slowness.

% . . Division in absolute slowness.
pmin=0.0399;
pmax=0.0800;
nbinp=40;
dpstep=(pmax-pmin)/nbinp;
pband=[(pmin:dpstep:pmax-dpstep);...
    (pmin+dpstep:dpstep:pmax)]';

% . . Division in back azimuth.
bzmin=0.0;
bzmax=360.0;
bzstep=5.0;
bzband=[(bzmin:bzstep:bzmax-bzstep);...
    (bzmin+bzstep:bzstep:bzmax)]';
nbinbz=size(bzband,1);

% . . Define slowness/baz grid.
npbz=nbinbz*nbinp;

bpbinmat=zeros(nfiles,nbinp,nbinbz);