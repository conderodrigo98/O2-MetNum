% Euler hacia atras

function [u,v,p,q] = EulerAtras (h,u0,v0,p0,q0,tol)
  
  pasos = 100/h;
  
	u = zeros(pasos+1,1) ;
	v = zeros(pasos+1,1) ;
	p = zeros(pasos+1,1) ;
	q = zeros(pasos+1,1) ;
	u(1) = u0 ;
	v(1) = v0 ;
	p(1) = p0 ;
	q(1) = q0 ;
  
  for i = 2:pasos+1  
    y = NR_Euler(u(i-1),v(i-1),p(i-1),q(i-1),h,tol);
    u(i) = y(1) ;
	  v(i) = y(2) ;
	  p(i) = y(3) ;
	  q(i) = y(4) ;
    
  end
  
end


