function [ fu, fv, fp, fq ] = evaluar(u,v,p,q)

fu = p ;
fv = q ;
fp = -2*v*p*q / (u^2+v^2+1) ;
fq = -2*u*p*q / (u^2+v^2+1) ;


end