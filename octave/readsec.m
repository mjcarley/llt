function [Re,al,cl,cd,iRe,sname]=readsec(file)

  dsize = 4 ;
  
  fid = fopen(file, "r") ;
  if ( fid < 0 )
    error(["cannot open file " file]) ;
  endif

  sname = fscanf(fid, "%[^\n]c") ;
  fscanf(fid, "%c", 1) ;

  nRe = fscanf(fid, "%d", 1) ;
  fscanf(fid, "%c", 1) ;
  iRe = fscanf(fid, "%d", nRe+1) ;
  fscanf(fid, "%c", 1) ;

  dat = fscanf(fid, "%f") ;

  Re = dat(1:dsize:end) ;
  al = dat(2:dsize:end) ;
  cl = dat(3:dsize:end) ;
  cd = dat(4:dsize:end) ;
  
  fclose(fid) ;
