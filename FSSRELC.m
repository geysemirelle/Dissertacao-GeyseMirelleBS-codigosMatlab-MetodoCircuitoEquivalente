%-----------------------------------------
% FSS ELC
%----------------------------------------------------
clear
clc
%
%--------------------------------------------------------------
% Entrada de dados codigo ELC
%--------------------------------------------------------------

f_inicial=0.5;
f_final=15;
n=0.01;
teta=0;
fi=0;
w1=0.9;
w2=0.9;
d1=2.35;
d2=13;
g1=1.8;
g2=1.8;
p=14.8;
er=4.3;
h=1.5;
%---------------------------------------------------------------

%Calculando a permissividade eletrica efetiva
x=10*h/p;
eff=er+(er-1)*(-1/((exp(x))^1.8));  


freq = f_inicial:n:f_final; %vetor de freq
n_f = length(freq); %iteracoes

lambda = 300./freq; %vetor de comprimento de onda lambda


XL1 = zeros(1,n_f); %pre-alocacao
XL2 = zeros(1,n_f);
BC1 = zeros(1,n_f);
BC2 = zeros(1,n_f);
y = zeros(1,n_f);
ct = zeros(1,n_f);
pr = zeros(1,n_f);
cr = zeros(1,n_f);

for ii = 1:n_f


 freq(ii) = f_inicial + (ii-1) * n;
 lambda(ii) = 300 / freq(ii);

X1=F(p,w1,lambda,fi); %calcula a reatancia indutiva X1
X2=F(p,w2,lambda,fi); %calcula a reatancia indutiva X2
X3=F(p,2*w2,lambda,fi); %calcula a reatancia indutiva X2

B1=4.*eff.*H(p,g1,lambda,teta); % calcula a susceptância da cruz
B2=4.*eff.*H(p,g2,lambda,teta); % calcula a susceptancia da cruz
BC1=0.37.*B1.*(d1./p);
BC2=4*0.75*(d2./p).*(B1.*B2./(B1+B2));
XL1=((d2-g1)/p)*X3; 
XL2=2.*(d2/p).*(X1.*X2./(X1+X2));

%reatancias
% (-1/BC1(ii))
% (-1/BC2(ii))
% XL1(ii)
% XL2(ii)



y(ii)=(1./(2.*XL1(ii)-(1./BC1(ii))))+(1./(XL2(ii)-(1./BC2(ii)))); %admitancia normalizada 

ct(ii)=10.*log10(4./(4+(y(ii).^2))); %coeficiente de transmissao 
end






%comparação 
% Frequência
f2= xlsread('quartacelula2.xlsx',1,'A2:A1002'); %frequencia teste 1
% Transmissão
x2= xlsread('quartacelula2.xlsx',1,'B2:B1002');


  %Gráficos
 figure(1)
 plot(freq,ct,f2,x2,'Linewidth',2); 
 hold on
 legend(['Análise numérica'],['Simulação computacional'])


% Considerando indutância TE incidente e Susceptância para TM incidente


%Calculando para Xte
%Calcula F(p,w,lambda,teta)

function Xte = F(p,w,lambda,teta)
A=((p.*(cos(teta)))./lambda);
B= log(csc((pi.*w)./(2.*p)));
Xte=A.*(B+G(p,w,lambda,teta));
end

%Calculando a funcao G, conforme referencia 
function funcaoG = G(p,w,lambda,teta)
beta=sin(pi.*w./(2.*p));
%Calculo de Cte1+ e Cte1-
D=(p.*sin(teta))./lambda;
E=(p.^2)./(lambda.^2);
Cte1_p=(1./(sqrt((D+1).^2-(E))))-1;
Cte1_n=(1./(sqrt((D-1).^2-(E))))-1;
F=(0.5.*((1-beta.^2).^2).*((1-(beta./2).^2).*(Cte1_p+Cte1_n)+4.*(beta.^2).*Cte1_p.*Cte1_n));
G=(1-(beta./2).^2)+(beta.^2).*(1+((beta.^2)./2)-((beta.^4)./8)).*(Cte1_p+Cte1_n)+2.*(beta.^6).*Cte1_p.*Cte1_n;
funcaoG=F./G;
end


%Calculando para Btm

function Btm = H(p,g,lambda,fi)
A1=((p.*(cos(fi)))./lambda);
B1= log(csc((pi.*g)./(2.*p)));
Btm=A1.*(B1+I(p,g,lambda,fi));
end

%Calculando a funcao G, conforme referencia 
function funcaoB = I(p,g,lambda,fi)
beta=sin(pi.*g./(2.*p));
%Calculo de Cte1+ e Cte1-
D1=(cos(fi));
E1=(p.^2)./(lambda.^2);
Ctm1_p=(1./(sqrt(1-(D1.*E1))))-1;
Ctm1_n=Ctm1_p;
F1=(0.5.*((1-beta.^2).^2).*((1-(beta./2).^2).*(Ctm1_p+Ctm1_n)+4.*(beta.^2).*Ctm1_p.*Ctm1_n));
G1=(1-(beta./2).^2)+(beta.^2).*(1+((beta.^2)./2)-((beta.^4)./8)).*(Ctm1_p+Ctm1_n)+2.*(beta.^6).*Ctm1_p.*Ctm1_n;
funcaoB=F1./G1;
end

%-----------------------------------------
