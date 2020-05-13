% 
  clear all;
  format compact;

%  dir='/Users/bostock/PGC/';
  dir='/fiore2_shared/shared/4Genevieve/PGC/';
  infile='list.COR';

  rlength=40;
%  ndecimate=10;
  ndecimate=2;
  hlength=200;
  
  ndt=0.025*ndecimate
  dtap=3.0;
  lf=0.025;
  hf=0.75;
  tshift=0.0;
  nt=round(tshift/ndt);
  nyq=0.5/(ndecimate*0.025);
  wn=[lf/nyq,hf/nyq];
  [b,a]=butter(2,wn);
  nbin=100
  difftimes;
  efac=180/(pi*6371);
