## generate a basic rectangular wing
ch = 0.3048 ;
b  = 30*ch ;

np = 60 ;

yb = linspace(-b/2,b/2,np+1)' ;
xb = [0*yb yb 0*yb] ;
##xb(:,1) += 0.12*ch ;

xw = xb  ;
xw(:,1) += 0.75*ch ;

xc = xb(1:end-1,:) + 0.5*diff(xb) ;
nc = xc*0 ; nc(:,3) = 1 ;
cc = xc*0 ; cc(:,1) = 1 ;
chc = ch*ones(size(xc(:,1))) ;

fid = fopen("wing.llt", "w") ;

fprintf(fid, "LINE %d\n", np) ;

dat = [xc nc cc chc] ;
fprintf(fid, "%e %e %e %e %e %e %e %e %e %e\n", dat') ;
dat = [xb xw] ;
fprintf(fid, "%e %e %e %e %e %e\n", dat') ;

fclose(fid) ;
