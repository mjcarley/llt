function G=elwing(a, b, U, al, y)

  ## Circulation on an elliptical wing with root chord a and span b at speed
  ## U and incidence al, -b/2 <= y <= b/2

  ## aspect ratio
  AR = 4*b/pi/a ;
  ## maximum circulation
  Gtm = 2*b*U^2/(1+AR/2)*al ;
  G = Gtm*sqrt(1-4*(y/b).^2) ;
  

