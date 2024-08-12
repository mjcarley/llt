function mksecfile(file, name, Re, al, cl, cd)

  ## Generate an LLT formatted section data file from plain text input
  [uRe,iRe] = unique(Re) ;

  fid = fopen(file, "w") ;

  fprintf(fid, "%s\n", name) ;
  fprintf(fid, "%d:", length(iRe)) ;
  fprintf(fid, " %d", iRe-1) ;
  fprintf(fid, " %d\n", length(Re)) ;

  dat = [Re al cl cd] ;

  fprintf(fid, "%e %e %e %e\n", dat') ;
  
  fclose(fid) ;
