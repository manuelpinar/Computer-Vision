%Generación media espacial de la corriente oscura (DCNU), varianza del DCNU y varianza temporal%
format long;
close all;
load(sprintf('10212-iproc%02da.mat',1));
M = size(MUDc,3);
EIMUDc = zeros(M,1);
VIMUDc = zeros(M,1);
EIS2Dc = zeros(M,1);
for k = 1:M
    % calculo de promedios y varianzas espaciales
    mu = MUDc(:,:,k);
    s2 = S2Dc(:,:,k);
    EIMUDc(k) = mean2(mu);
    VIMUDc(k) = var(mu(:));
    EIS2Dc(k) = mean2(s2);
end

elec=menu('TIPO DE AJUSTE','Sin Outliers', 'Con Outliers')

% Inyección de OUTLIERS
if elec==2
    load(sprintf('10212-outliers%02d.mat',1));
    EIMUDc = EIMUDc + outl1;
    VIMUDc = VIMUDc + outl2;
    EIS2Dc = EIS2Dc + outl3;
end

%Utilización de Polyfit para ajustar la parábola y la recta%
p=polyfit(EIMUDc,VIMUDc,2);

x=0:0.1:255;subplot(2,2,1);plot(x,p(1)*x.^2+p(2)*x+p(3));title('POLYFIT PARÁBOLA');
hold on;
plot(EIMUDc,VIMUDc,'+r');

p2=polyfit(EIMUDc,EIS2Dc,1);

[VIK, Ac, Var_ruido_aleatorio2, EDCNU, VDCNU]=parametros2(p(1),p(2),p(3),p2(1),p2(2));
fid = fopen('polyfit.txt', 'wt');
B = [VIK Ac Var_ruido_aleatorio2 EDCNU VDCNU]
fprintf(fid, 'VIK es %8.10f\n',B(1));
fprintf(fid, 'Ac es %8.10f\n',B(2));
fprintf(fid, 'Var_ruido_aleatorio2 es %8.10f\n',B(3));
fprintf(fid, 'EDCNU es %8.10f\n',B(4));
fprintf(fid, 'VDCNU es %8.10f\n',B(5));
fclose(fid)

x=0:0.1:255;subplot(2,2,2);plot(x,p2(1)*x+p2(2)),title('POLYFIT RECTA');
hold on;
plot(EIMUDc,EIS2Dc,'+r');

%B1%
%error con polinomio de orden 2%
Yestimadap=polyval(p,EIMUDc); 
Err_ajustep=abs(VIMUDc-Yestimadap);
Media_errorp=mean(Err_ajustep);
Max_errorp=max(Err_ajustep);
%error con polinomio de orden 1%
Yestimadap2=polyval(p2,EIMUDc);
Err_ajustep2=abs(EIS2Dc-Yestimadap2);
Media_errorp2=mean(Err_ajustep2);
Max_errorp2=max(Err_ajustep2);

%B2 Utilización de Polytool para ajustar la parábola y la recta
%Para el polinomio de orden 2%
polytool(EIMUDc,VIMUDc,2); %intervalos de confianza de a,b y c de P y errores de Yestimada respecto de Y real%
pause                     
residualsp=evalin('base','residualsp');
Int_confp=evalin('base','Int_confp');
par=evalin('base','par');
residualsp=abs(residualsp);
media_residualp=mean(residualsp);
max_residualsp=max(residualsp);
%Para el polinomio de orden 1%

polytool(EIMUDc,EIS2Dc,1); %intervalos de confianza de a,b de P2 y errores de Yestimada respecto de Y real%
pause                          
residualsp2=evalin('base','residualsp2');
Int_confp2=evalin('base','Int_confp2');
par2=evalin('base','par2');
residualsp2=abs(residualsp2);
media_residualp2=mean(residualsp2);
max_residualsp2=max(residualsp2);

[VIK, Ac, Var_ruido_aleatorio2, EDCNU, VDCNU]=parametros2(par(1),par(2),par(3),par2(1),par2(2));
fid = fopen('polytool.txt', 'wt');
B = [VIK Ac Var_ruido_aleatorio2 EDCNU VDCNU]
fprintf(fid, 'VIK es %8.10f\n',B(1));
fprintf(fid, 'Ac es %8.10f\n',B(2));
fprintf(fid, 'Var_ruido_aleatorio2 es %8.10f\n',B(3));
fprintf(fid, 'EDCNU es %8.10f\n',B(4));
fprintf(fid, 'VDCNU es %8.10f\n',B(5));
fclose(fid)

%B3 Utilización de RANSAC para ajustar la parábola y la recta%
 
%Ransac para ajuste de parábola
data(:,1)=EIMUDc;
data(:,2)=VIMUDc;
coef=[1 1 1];
k=500;
ransac=RANSACpolyfit(data,coef,k,0.6,0.8,1);
x=0:0.1:255;subplot(2,2,3);plot(x,ransac(1)*x.^2+ransac(2)*x+ransac(3)),title('RANSAC "PARÁBOLA"');
hold on;
plot(EIMUDc,VIMUDc,'+r');
%Ransac para ajuste de recta
data2(:,1)=EIMUDc;
data2(:,2)=EIS2Dc;
coef2=[1 1];
ransac2=RANSACpolyfit(data2,coef2,k,0.6,0.8,1);
x=0:0.1:255;subplot(2,2,4);plot(x,ransac2(1)*x+ransac2(2)),title('RANSAC "RECTA"');
hold on;
plot(EIMUDc,EIS2Dc,'+r');

[VIK, Ac, Var_ruido_aleatorio2, EDCNU, VDCNU]=parametros2(ransac(1),ransac(2),ransac(3),ransac2(1),ransac2(2));
fid = fopen('ransac.txt', 'wt');
B = [VIK Ac Var_ruido_aleatorio2 EDCNU VDCNU]
fprintf(fid, 'VIK es %8.10f\n',B(1));
fprintf(fid, 'Ac es %8.10f\n',B(2));
fprintf(fid, 'Var_ruido_aleatorio2 es %8.10f\n',B(3));
fprintf(fid, 'EDCNU es %8.10f\n',B(4));
fprintf(fid, 'VDCNU es %8.10f\n',B(5));
fclose(fid)

%error con ransac de orden 2%
Yestimadaransac=polyval(ransac,EIMUDc); 
Err_ajusteransac=abs(VIMUDc-Yestimadaransac);
Media_errorransac=mean(Err_ajusteransac);
Max_errorransac=max(Err_ajusteransac);

%error con ransac de orden 1%
Yestimadaransac2=polyval(ransac2,EIMUDc);
Err_ajusteransac2=abs(EIS2Dc-Yestimadaransac2);
Media_errorransac2=mean(Err_ajusteransac2);
Max_errorransac2=max(Err_ajusteransac2);

%B4
%LSQNONNEG DE ORDEN 2
A=[EIMUDc.^2 -2*EIMUDc ones(15,1)];
noneg=lsqnonneg(A,VIMUDc);
Yest=polyval(noneg,EIMUDc); 
Err_ajustenoneg=abs(VIMUDc-Yest);
Media_errornoneg=mean(Err_ajustenoneg);
Max_errornoneg=max(Err_ajustenoneg);
x=0:0.1:255;figure;subplot(2,2,1);plot(x,noneg(1)*x.^2-2*noneg(2)*x+noneg(3)),title('LSQNONNEG PARÁBOLA');
hold on;
plot(EIMUDc,VIMUDc,'+r');

%LSQNONNEG DE ORDEN 1
A2=ones(15,2);
A2(:,1)=EIMUDc;
noneg2=lsqnonneg(A2,EIS2Dc);
Yest2=polyval(noneg2,EIMUDc); 
Err_ajustenoneg2=abs(EIS2Dc-Yest2);
Media_errornoneg2=mean(Err_ajustenoneg2);
Max_errornoneg2=max(Err_ajustenoneg2);
x=0:0.1:255;subplot(2,2,2);plot(x,noneg2(1)*x+noneg2(2)),title('LSQNONNEG RECTA');
hold on;
plot(EIMUDc,EIS2Dc,'+r');

[VIK, Ac, Var_ruido_aleatorio2, EDCNU, VDCNU]=parametros(noneg(1),noneg(2),noneg(3),noneg2(1),noneg2(2));
fid = fopen('lsqnonneg.txt', 'wt');
B = [VIK Ac Var_ruido_aleatorio2 EDCNU VDCNU]
fprintf(fid, 'VIK es %8.10f\n',B(1));
fprintf(fid, 'Ac es %8.10f\n',B(2));
fprintf(fid, 'Var_ruido_aleatorio2 es %8.10f\n',B(3));
fprintf(fid, 'EDCNU es %8.10f\n',B(4));
fprintf(fid, 'VDCNU es %8.10f\n',B(5));
fclose(fid)

%B5
%Ransac manipulado para ajuste de parábola
data(:,1)=EIMUDc;
data(:,2)=VIMUDc;
coef=[1 -2 1];
ransac_pos=RANSACpolyfit(data,coef,k,0.5,0.7,0);
x=0:0.1:255;subplot(2,2,3);plot(x,ransac_pos(1)*x.^2+ransac_pos(2)*x+ransac_pos(3)),title('RANSAC PARÁBOLA MANIPULADO');
hold on;
plot(EIMUDc,VIMUDc,'+r');

%Ransac manipulado para ajuste de recta
data2(:,1)=EIMUDc;
data2(:,2)=EIS2Dc;
ransac2_pos=RANSACpolyfit(data2,coef2,k,0.5,0.7,0);
x=0:0.1:255;subplot(2,2,4);plot(x,ransac2_pos(1)*x+ransac2_pos(2)),title('RANSAC RECTA MANIPULADO');
hold on;
plot(EIMUDc,EIS2Dc,'+r');

[VIK, Ac, Var_ruido_aleatorio2, EDCNU, VDCNU]=parametros(ransac_pos(1),ransac_pos(2),ransac_pos(3),ransac2_pos(1),ransac2_pos(2));
fid = fopen('ransacpos.txt', 'wt');
%Z = ['VIK','Ac','Var_ruido_aleatorio2','EDCNU','VDCNU'];
B = [VIK Ac Var_ruido_aleatorio2 EDCNU VDCNU];
fprintf(fid, 'VIK es %8.10f\n',B(1));
fprintf(fid, 'Ac es %8.10f\n',B(2));
fprintf(fid, 'Var_ruido_aleatorio2 es %8.10f\n',B(3));
fprintf(fid, 'EDCNU es %8.10f\n',B(4));
fprintf(fid, 'VDCNU es %8.10f\n',B(5));
fclose(fid)

%error con ransac manipulado de orden 2%
Yestimadaransac_pos=polyval(ransac_pos,EIMUDc); 
Err_ajusteransac_pos=abs(VIMUDc-Yestimadaransac_pos);
Media_errorransac_pos=mean(Err_ajusteransac_pos);
Max_errorransac_pos=max(Err_ajusteransac_pos);

%error con ransac manipulado de orden 1%
Yestimadaransac2_pos=polyval(ransac2_pos,EIMUDc);
Err_ajusteransac2_pos=abs(EIS2Dc-Yestimadaransac2_pos);
Media_errorransac2_pos=mean(Err_ajusteransac2_pos);
Max_errorransac2_pos=max(Err_ajusteransac2_pos);

% 3 VALORES DE INCERTIDUMBRE

% Apartado A
%Utilizamos los coeficientes obtenidos en el ransac manipulado por ser los
%que ofrecen una media de error más baja
 
varK=par(1);
Ei_dcAc=(-par(2))/(2*par(1));
vardcAc=(par(3)-(par(2)^2/(4*par(1))));
Ac_pol=par2(1);
sigma=par2(2);
t=3;
Dc=0:1:255;

for i=1:256
    alfa=(1-(t^2)*varK);
    beta=(((t^2)*Ac_pol)+2*(Dc(i)-Ei_dcAc));
    gamma=(Dc(i)-Ei_dcAc)^2-(t^2)*((vardcAc)^2+Ei_dcAc*(Ac_pol)+sigma);
    pol=[alfa -beta gamma];
    polinomio=roots(pol);
    polinomio=(polinomio)';
    intensidades(i,1)=polinomio(2); %Creo una matriz con los valores de intensidades Ia y Ib
    intensidades(i,2)=polinomio(1);
    dig(i,1)=polinomio(2)+Ei_dcAc; 
    dig(i,2)=polinomio(1)+Ei_dcAc;
end
figure;subplot(2,2,1);plot(Dc,dig(:,1),'-r'),title('Intervalo de confianza exacto para t=3');
hold on 
plot(Dc,dig(:,2),'-b');

resultadoA=[[100:110]' dig(100:110,1:2)]; %Matriz de resultados Niveles de incertidumbre para valores de intensidad entre 100 y 110

%APARTADO B
% INTERVALO DE CONFIANZA APROXIMADO PARA t=1 Y t=3

for i=1:256
    sigmaint=sqrt((varK*(Dc(i)-Ei_dcAc)^2)+(vardcAc)+((Dc(i)-Ei_dcAc)*Ac_pol)+(Ei_dcAc*Ac_pol)+sigma);
    dig_aproxT1(i,1:2)=[Dc(i)-sigmaint Dc(i)+sigmaint];
    dig_aproxT3(i,1:2)=[Dc(i)-t*sigmaint Dc(i)+t*sigmaint];
end
subplot(2,2,2);plot(Dc,dig_aproxT1(:,1),'-r'),title('Intervalo de confianza aproximado para t=1');
hold on 
plot(Dc,dig_aproxT1(:,2),'-b');
subplot(2,2,3);plot(Dc,dig_aproxT3(:,1),'-r'),title('Intervalo de confianza aproximado para t=3');
hold on 
plot(Dc,dig_aproxT3(:,2),'-b');
resultadoB1=[[100:110]' dig_aproxT1(100:110,1:2)]; %Matriz de resultados Niveles de incertidumbre para valores de intensidad entre 100 y 110 con t=1
resultadoB3=[[100:110]' dig_aproxT3(100:110,1:2)]; %Matriz de resultados Niveles de incertidumbre para valores de intensidad entre 100 y 110 con t=3

%APARTADO C
%DIFERENCIA ENTRE INTERVALO EXACTO Y APROXIMADO CON t=3
dif_intervalos=abs(dig-dig_aproxT3);  
subplot(2,2,4);plot(Dc,dif_intervalos(:,1),'-r'),title('Diferencia de intervalos de confianza aproximado para t=3');
hold on 
plot(Dc,dif_intervalos(:,2),'-b');
resultadoC=[[100:110]' dif_intervalos(100:110,1:2)]; %Matriz de resultados de diferncia  Niveles de incertidumbre para valores de intensidad entre 100 y 110 con t=3
