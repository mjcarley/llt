function dat=readprop(file)

  fid = fopen(file, "r") ;

  ## read and ignore first three lines
  fscanf(fid, "%[^\n]s") ; fscanf(fid, "%c", 1) ;
  fscanf(fid, "%[^\n]s") ; fscanf(fid, "%c", 1) ;
  fscanf(fid, "%[^\n]s") ; fscanf(fid, "%c", 1) ;

  dat = [] ;
  [buf,n] = fscanf(fid, "%[^\n]s") ; fscanf(fid, "%c", 1) ;
  ii = find(buf == "\t") ;
  buf = buf(ii(1):end) ;
  lline = sscanf(buf, "%f") ;
  dat = [dat; lline'] ;

  while ( lline(1) ~= 1 )
    [buf,n] = fscanf(fid, "%[^\n]s") ; fscanf(fid, "%c", 1) ;
    ii = find(buf == "\t") ;
    buf = buf(ii(1):end) ;
    lline = sscanf(buf, "%f") ;
    dat = [dat; lline'] ;
  endwhile
  
  fclose(fid) ;

  ## data columns on output are:
  ##
  ##  1: r/R
  ##  2: c/R
  ##  3: ph (pitch)
  ##  4: ph (pitch)
  ##  5: ph (pitch)?
  ##  6: sw (sweep)
  ##  7: th (thickness)
  ##  8: tw (twist)
  ##  9: thmax (maximum thickness)
  ## 10: xsec (cross-section)
  ## 11: zhigh
  ## 12: cgy
  ## 13: cgz
  
