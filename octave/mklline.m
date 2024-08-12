## input file containing propeller blade data
fname = "APC10x5E.txt" ;
dat = readprop(fname) ;
## 5in radius
R = 5*25.4/1000 ;

## generate lifting line representation
r  = dat(:,1)*R ;
ch = dat(:,2)*R ;
ph = dat(:,3) ;
sw = dat(:,6) ;
th = dat(:,7) ;
tw = dat(:,8) ;

tw *= pi/180 ;
## output control points x n c ch
## output bound points x w
rc = r(1:end-1) + diff(r)*0.5 ;
chc = interp1(r, ch, rc, "cubic") ;
tc = interp1(r, tw, rc, "cubic") ;

xc = [0*tc rc 0*tc] ;
#nc = [cos(tc) 0*tc -sin(tc)] ;
##cc = [sin(tc) 0*tc  cos(tc)] ;
##c = [sin(tw) 0*tw  cos(tw)] ;
cc = [cos(tc) 0*tc -sin(tc)] ;
nc = [sin(tc) 0*tc  cos(tc)] ;
c = [cos(tw) 0*tw  -sin(tw)] ;

fid = fopen("apc-prop.llt", "w") ;
fprintf(fid, "LINE %d\n", length(rc)) ;

dat = [xc nc cc chc] ;
fprintf(fid, "%e %e %e %e %e %e %e %e %e %e\n", dat') ;

xb = [0*r    r 0*r] ;
xw = xb + [0.75*ch.*c(:,1) 0.75*ch.*c(:,2) 0.75*ch.*c(:,3)] ; 

dat = [xb xw] ;
fprintf(fid, "%e %e %e %e %e %e\n", dat') ;

fclose(fid) ;
