%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% FBxLMS Distribuido Colaborativo INCREMENTAL SIN RESTRICCIONES
clear 
clc
disp('');
disp('FBPFxLMS_Colaborativo Incremental');
rootfolder = 'E:\Doctorado\Doctorado TELECO\Tesis\Pruebas MATLAB\Mis Pruebas';
addpath([rootfolder,'\Datos\Señales']);
addpath([rootfolder,'\Datos\Caminos']);
%% Configuración de la red
B=256;
ITE=100000;    % nº iteraciones 
L=150;
tam_fft=2*B;   % tamaño fft

%% señal de entrada
load x3;
in=1.5*x3(1:ITE,1);
x=in-ones(ITE,1)*mean(in);
N=length(x); 

%% Cargar Sala 
load ir_map; 
sala=ir_map/(max(max(max(max(ir_map)))));
fs=2000;
%% Definición de los datos 
%% RED 4 NODOS 
alt=[20 21 22 23];
mic=[61 60 59 58]; 
mu=0.002; 
beta=0.001;
%% RED 10 NODOS
% alt=[15,16,17,18,19,20,21,22,23,24];
% mic=[18,17,16,15,14,13,12,11,10,9];
% mu=0.0006;% 
% beta=0.001;
%% RED 16 NODOS
% alt=[ 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
% mic=[71 70 69 68 67 66 65 64 63 62 61 60 59 58 57 56];
% alt=[15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30];
% mic=[18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3];
% mu=0.001;
% beta=0.000;



%% Fuente de ruido
ref=69;
%% Número de fuentes de ruido, altavoces y micros
I=length(ref);
J=length(alt);
K=length(mic);

%% Definición Nodos
num_nodos=K;
nodos=[alt;
       mic];
mu=mu/num_nodos; 
   
%% Cálculo Caminos Acústicos  
% Para los nodos
for n=1:num_nodos
    c=sala(alt(1:end),:,nodos(2,n))';
    csec(:,:,n)=c;
end
% Para el Sistema acústico
M=length(sala(1,:,1,1));   
C_aux=sala(ref(1:I),:,mic(1:K));
C=permute(C_aux,[2 3 1]);
if I==1
   pri=C;
else
   pri=reshape(C,M,I*J);
end 
C_aux=sala(alt(1:J),:,mic(1:K));
C=permute(C_aux,[2 3 1]);
CSEC=reshape(C,M,J*K);

% Cálculo de señal deseada
d=zeros(N,K);
for k=1:K
    d(:,k)=filter(pri(:,k),1,x);
end

e=zeros(length(d),num_nodos);
IN_e=zeros(B,num_nodos);
num=floor(N/B);

%% Configuración del algoritmo
% Particiones
F=ceil(M/B);
P=ceil(L/B);
% Caminos en frecuencia
c_aux=[csec;zeros(B*F-M,num_nodos,num_nodos)];
cp=reshape(c_aux,B,F,num_nodos,num_nodos);
C=fft(cp,tam_fft,1);

C_aux=[CSEC;zeros(B*F-M,J*K)];
Cp=reshape(C_aux,B,F,J*K);
CSEC_F=fft(Cp,tam_fft,1);

% Inicialización de variables
W=zeros(num_nodos*tam_fft,P,num_nodos);
V=zeros(num_nodos*tam_fft,P,num_nodos);
buffx=zeros(tam_fft,1);
x_f=zeros(tam_fft,F,num_nodos);
x_p=zeros(tam_fft,P,num_nodos);

IN_y=zeros(B,num_nodos);
buffy=zeros(tam_fft,num_nodos);
y_f=zeros(tam_fft,F,num_nodos);

control=zeros(length(d),num_nodos);

O_est=zeros(num_nodos*tam_fft,P);

%% Ejemplos de uso de beta, ymax, alfa:
beta=ones(num_nodos,1)*beta;
vectorbeta=ones(L,1)*beta(:)';
vectorbeta=vectorbeta(:);

ymax=ones(J,1)*1;
Ymax=repmat(ymax',B,1);



alfa=ones(J,1)*1;

%% Comienzo de la simulación %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic 
for num_b=1:num
    
  % Señal entrada
  IN=x(B*num_b-B+1:B*num_b);
  % Buffer de tamaño bloque de la señal de entrada   
  buffx=[buffx(B+1:tam_fft);IN];
    
  for nodo=1:num_nodos 
    
    % Para simular red incremental
    if nodo==1
        W_ant=W(:,:,num_nodos);
    else
        W_ant=W(:,:,nodo-1);
    end  
    
    % Utilizo sólo los coeficientes del nodo para el filtrado adaptativo
    W_nodo=W(1+tam_fft*(nodo-1):tam_fft*nodo,:,nodo);
    
    % Filtrado adaptativo 
    x_p(:,2:end,nodo)=x_p(:,1:P-1,nodo);
    x_p(:,1,nodo)=buffx;  
    X_p=fft(x_p,tam_fft);
    Y=X_p(:,:,nodo).*W_nodo;                  
    y=ifft(Y); 
    yb=sum(y,2);     
    IN_y(:,nodo)=yb(B+1:tam_fft,:); 
    
    E=fft([zeros(B,1);IN_e(:,nodo)],tam_fft); 
    Et=repmat(E,num_nodos,P);      
    
    phi=Et.*conj(V(:,:,nodo));
          
    for n=1:num_nodos
       step=ifft(phi(tam_fft*(n-1)+1:tam_fft*n,:));
       O=fft([step(1:B,:);zeros(B,P)],tam_fft);
       O_est(tam_fft*(n-1)+1:tam_fft*n,:)=O;
    end
    % Actualización de coeficientes
    W(:,:,nodo)=W_ant-2*mu.*(((beta(nodo)/num_nodos)*W(:,:,nodo))+O_est);
    
    % Filtrado estima
    x_f(:,2:end,nodo)=x_f(:,1:F-1,nodo);
    x_f(:,1,nodo)=buffx;  
    X_f=fft(x_f,tam_fft);
    for j=1:J
       V_aux=X_f(:,:,nodo).*C(:,:,j,nodo);
       Vaux=sum(V_aux,2);
       V_a(1+tam_fft*(j-1):tam_fft*j,:)=Vaux;
    end
    V(:,2:end,nodo)=V(:,1:P-1,nodo);
    V(:,1,nodo)=V_a;    
    
    %% Segunda vuelta para el Reescalado
        if (ymax(nodo)~=0)       
          k=find(abs(IN_y(:,nodo))>ymax(nodo));
          if k
            IN_y_ant=IN_y; 
            IN_y(k,nodo)=IN_y(k,nodo)*alfa(nodo).*(ymax(nodo)./abs(IN_y(k,nodo)));
            W(1+tam_fft*(nodo-1):tam_fft*nodo,:,nodo)=W(1+tam_fft*(nodo-1):tam_fft*nodo,:,nodo)*alfa(nodo)*mean(ymax(nodo)./abs(IN_y_ant(k,nodo)));
          end
        end         
  end
  
  % Difusión de los coeficientes del último nodo
%   W=repmat(W(:,:,end),1,1,num_nodos);
   
    
 %% SISTEMA ACÚSTICO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  buffy=[buffy(B+1:2*B,:);IN_y];
  y_f(:,2:end,:)=y_f(:,1:F-1,:);
  y_f(:,1,:)=buffy;  
  YF=fft(y_f);

  Yf=zeros(tam_fft,F,K);
  for k=1:K
    for j=1:J
       Yf(:,:,k)=Yf(:,:,k)+YF(:,:,j).*CSEC_F(:,:,k+K*(j-1));
    end
  end
    
  yf=ifft(Yf); 
  yf=sum(yf,2);    
  OUT=yf(B+1:end,:);
       
  % Señal deseada
  db=d(B*num_b-B+1:B*num_b,:); 
  IN_e=db+OUT; 
  e(B*num_b-B+1:B*num_b,:)=IN_e; 
  control(B*num_b-B+1:B*num_b,:)=IN_y;
 
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end 
toc 
%% Fin simulación %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% GRAFICAS
dpot=calpot(d,0.9995);
epot=calpot(e,0.9995);

ATEp=epot./dpot;
ATEp_db=10*log10(ATEp);
t=1/fs*(1:length(ATEp))';
% figure;
hold all;
plot(t,ATEp_db(:,2));
axis([0,t(end),-20,5]);
ylabel('Noise Reduction (dB)')
xlabel('Time (s)')

