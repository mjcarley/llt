function [x,y,z,G]=readwake(file)

  psize = 5 ;
  fid = fopen(file, "r") ;

  dat = fscanf(fid, "%d", 2) ;

  nt = dat(1) ; ne = dat(2) ;

  dat = fscanf(fid, "%f", ne*nt*psize) ;

  x = reshape(dat(1:psize:end), ne, nt)' ;
  y = reshape(dat(2:psize:end), ne, nt)' ;
  z = reshape(dat(3:psize:end), ne, nt)' ;
  G = reshape(dat(4:psize:end), ne, nt)' ;
  
  fclose(fid) ;
