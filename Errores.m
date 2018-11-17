% Errores

clear all, close all, clc

u0 = 1 ;
v0 = 1 ; 
p0 = 0.01 ; 
q0 = 0.01 ;
tol = 1e-6 ;


hini = 0.01 ;
hobj = 1 ;

vectorH = [hini : 0.01 : hobj] ;

EmaxR = zeros(length(vectorH),1) ;
EmaxH = zeros(length(vectorH),1) ;
EmaxBE = zeros(length(vectorH),1) ;
EmaxEA = zeros(length(vectorH),1) ;
EmaxT = zeros(length(vectorH),1) ;




for i = 1:length(vectorH)

  h = vectorH(i) ;
  [uR,vR,pR,qR] = Runge(h,u0,v0,p0,q0,tol);
  [uH,vH,pH,qH] = Heun(h,u0,v0,p0,q0);
  [uBE,vBE,pBE,qBE] = EulerAtras (h,u0,v0,p0,q0,tol);
  [uFE,vFE,pFE,qFE] = EulerAdelante (h,u0,v0,p0,q0) ;
  [uT,vT,pT,qT] = Trapecio (h,u0,v0,p0,q0,tol) ;
  
  pasos = 100/h;  
  t = linspace(0,100,pasos+1);
  y = lsode("f",[ 1; 1; 0.01; 0.01],t);
  
  UV_R = uR.*vR;
  UV_H = uH.*vH;
  UV_BE = uBE.*vBE;
  UV_FE = uFE.*vFE;
  UV_T = uT.*vT;
  
  UVlsode = y(1:pasos+1,1).*y(1:pasos+1,2);
  
  ErroresR = abs(UVlsode-UV_R) ;
  ErroresH = abs(UVlsode-UV_H) ;
  ErroresBE = abs(UVlsode-UV_BE) ;
  ErroresFE = abs(UVlsode-UV_FE) ;
  ErroresT = abs(UVlsode-UV_T) ;
  
  
  EmaxR(i) = max(ErroresR) ;
  EmaxH(i) = max(ErroresH) ;
  EmaxBE(i) = max(ErroresBE) ;
  EmaxFE(i) = max(ErroresFE) ;
  EmaxT(i) = max(ErroresT) ;
  
end

plotfontsize = 22; 
lw = 1.2;
ms = 5.5;

plot(vectorH,EmaxR,'b-','linewidth',lw,'markersize',ms)
tit = title('Estudio error maximo - h optimo') ;
set(tit, "FontSize",plotfontsize) ;
labx = xlabel('Paso h') ;
laby = ylabel('Error maximo') ;
set(labx, "FontSize",plotfontsize) ; set(laby, "FontSize",plotfontsize) ;
 hold on
plot(vectorH,EmaxH, 'r' ,'linewidth',lw)
plot(vectorH,EmaxBE, ' k' ,'linewidth',lw)
plot(vectorH,EmaxFE, ' g','linewidth',lw )
plot(vectorH,EmaxT, ' m','linewidth',lw )

legend(' Runge' ,' Heun' , ' BE' , ' FE' ,' T' )

print(['errores'], ' -dpng' ) ;

fileDatos = fopen('./latex.tex' , 'w');
fprintf(fileDatos, ['\\begin{table}[H] \n' ] )
fprintf(fileDatos, ['\\centering \n' ] )
fprintf(fileDatos, ['\\begin{tabular}{cccccc} \n' ] )
fprintf(fileDatos, ['h & EmaxR & EmaxH & EmaxBE & EmaxFE & EmaxT \\\\ \\toprule \n' ] )
for i = 1:10

fprintf(fileDatos,[' %i & %3i & %3i & %3i & %3i & %3i \\\\ \\midrule \n' ], [ vectorH(i), EmaxR(i), EmaxH(i), EmaxBE(i), EmaxFE(i), EmaxT(i) ] ) 

end
fprintf(fileDatos, ['\\end{tabular} \n' ] )
fprintf(fileDatos, ['\\end{table} \n' ] )

fclose(fileDatos);



