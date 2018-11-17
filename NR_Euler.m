function answ = NR_Euler(un,vn,pn,qn,h,tol)
  
  function J = Jacobiano(u,v,p,q)
        J = zeros(4,4);
        J(1,1) = 1;
        J(1,2) = 0;
        J(1,3) = -h;
        J(1,4) = 0;
        J(2,1) = 0;
        J(2,2) = 1;
        J(2,3) = 0;
        J(2,4) = -h;
        J(3,1) = -h*4*u*v*p*q/((u^2+v^2+1)^2);
        J(3,2) = 2*h*p*q*(-v^2+u^2+1)/((u^2+v^2+1)^2);
        J(3,3) = 1+2*h*v*q/(u^2+v^2+1);
        J(3,4) = 2*h*v*p/(u^2+v^2+1);
        J(4,1) = 2*h*p*q*(-u^2+v^2+1)/((u^2+v^2+1)^2);
        J(4,2) = -h*4*u*v*p*q/((u^2+v^2+1)^2);
        J(4,3) = 2*h*u*q/(u^2+v^2+1);
        J(4,4) = 1+2*h*u*p/(u^2+v^2+1);       
  end
  
  anterior = [un;vn;pn;qn];
  actual = [un;vn;pn;qn];
  k=0;
  
  while tol > norm(anterior-actual) && k < 1000000
    anterior = actual;
    
    d = Jacobiano(anterior(1),anterior(2),anterior(3),anterior(4))\-evaluarG(anterior(1),anterior(2),anterior(3),anterior(4));
    actual = anterior + d;
    
    k++;
  end
  
  answ = actual;
  
  function res = evaluarG(u,v,p,q)
    res = zeros(4,1);
    res(1) = u - un - h*p;
    res(2) = v - vn -h*q;
    res(3) = p - pn + 2*h*(v*p*q)/(u^2+v\2+1);
    res(4) = q - qn + 2*h*(u*p*q)/(u^2+v\2+1);
  end
  
end