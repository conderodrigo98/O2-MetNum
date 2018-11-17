% Runge

function [u,v,p,q] = Runge(h, u0, v0, p0, q0)

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
   
    [fu,fv,fp,fq] = evaluar( u(i-1), v(i-1), p(i-1), q(i-1) ) ; 
    K1 = h * [ fu , fv , fp , fq ] ;
    
    [fu,fv,fp,fq] = evaluar( u(i-1) + h/2*K1(1), v(i-1) + h/2*K1(2), p(i-1) + h/2*K1(3), q(i-1) + h/2*K1(4) ) ;
    K2 = h * [ fu , fv , fp , fq ] ;
    
    [fu,fv,fp,fq] = evaluar( u(i-1) + h/2*K2(1), v(i-1) + h/2*K2(2), p(i-1) + h/2*K2(3), q(i-1) + h/2*K2(4) ) ;
    K3 = h * [ fu , fv , fp , fq ] ;
   
    [fu,fv,fp,fq] = evaluar( u(i-1) + h*K3(1), v(i-1) + h*K3(2), p(i-1) + h*K3(3), q(i-1) + h*K3(4) ) ;
    K4 = h * [ fu , fv , fp , fq ] ;
    
    u(i) = u(i-1) + 1/6 * ( K1(1) + 2*K2(1) + 2*K3(1) + K4(1) ) ;
    v(i) = v(i-1) + 1/6 * ( K1(2) + 2*K2(2) + 2*K3(2) + K4(2) ) ;
    p(i) = p(i-1) + 1/6 * ( K1(3) + 2*K2(3) + 2*K3(3) + K4(3) ) ;
    q(i) = q(i-1) + 1/6 * ( K1(4) + 2*K2(4) + 2*K3(4) + K4(4) ) ;
  
  end
  
  
end