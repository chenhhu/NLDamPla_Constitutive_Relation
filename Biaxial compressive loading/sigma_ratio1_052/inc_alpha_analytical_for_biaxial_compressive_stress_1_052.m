% This code is developed to give the analytical solution of the 
% nonlocal continuum damage-plasticity model, see Chen et al., 2024 (Journal). 
% Stress state: Uniaxial compression stree state
% test data-Kupfer and Gerstle, 1973
clear all
close all
clc
%%
% csvread(filename,R1,C1,[R1 C1 R2 C2]) extracted from experimental datas 
y1=csvread('Ratio1_52.csv',0,0,[0,0,12,1])';
xdata=y1(1,:);
ydata=y1(2,:);
%%
syms chi_p E c_f alpha_p c_0 h_p beta kappa kappa_i sigma1 A B D
% local equivalent strain
kappa_fun0=sqrt( 2*chi_p/E*( c_f*alpha_p+(c_f-c_0)/h_p*(exp(-h_p*alpha_p)-1) ) );
kappa_fun1=str2func(['@(chi_p,E,c_f,alpha_p,c_0,h_p)',vectorize(kappa_fun0)]); % Construct function handle from function 
% damage evolution
D_fun0=1-exp(-beta*(kappa-kappa_i));
D_fun1=str2func(['@(beta,kappa,kappa_i)',vectorize(D_fun0)]); % Construct function handle from function 
% The following is the yield function under Uniaxial compression
f_fun0 = sqrt(1.5008*sigma1^2/6)+A/3*sigma1*1.52-B*(1-chi_p*D)*(c_f-(c_f-c_0)*exp(-h_p*alpha_p))==0; 
sigma1_fun0=solve(f_fun0,sigma1); % sigma1 = f(alpha_p)
sigma1_fun1=str2func(['@(A,B,chi_p,D,c_f,c_0,h_p,alpha_p)',vectorize(sigma1_fun0(1))]); % Construct function handle from function "d2sigma1".
%% Initiation of material parameters. 
i=0;
E=31000; %[MPa]
v=0.2;
c_0=4.0569;
theta=pi/180*40;
A=6*sin(theta)/sqrt(3)/(3+sin(theta));
B=6*cos(theta)/sqrt(3)/(3+sin(theta));
C=6*sin(theta)/sqrt(3)/(3+sin(theta));

par=[1.0000   12.7673   58.8665   88.0123  478.3554];
chi_p=par(1);
sigma_c=par(2);                        
kappa_i=sigma_c/E;
c_f=par(3);
h_p=par(4);
beta=par(5);
                                          
%% Load by increasing d
for alpha_p=0.000001:0.0005:0.075
    i=i+1;
    kappa=kappa_fun1(chi_p,E,c_f,alpha_p,c_0,h_p);
    if kappa<kappa_i
        D(i)=0;
    else
        D(i)=D_fun1(beta,kappa,kappa_i);
    end
    sigma1(i)=sigma1_fun1(A,B,chi_p,D(i),c_f,c_0,h_p,alpha_p); % Axial stress
    N=-1.48/sqrt(9.0048)+C/3; % Axial plastic flow direction
    strain1(i)=alpha_p/B*N+sigma1(i)/E-0.52*v*sigma1(i)/E; % Axial strain 
    ana_alpha_p(i)=alpha_p;
    ana_sigma1(i)=(3*B/(1.52*A-sqrt(2.2512)))*(1-chi_p*D(i))*(c_f-(c_f-c_0)*exp(-h_p*alpha_p));     
end  
%% 
figure(1)
plot(xdata,ydata,'ko',[0,strain1],[0,sigma1],'b-');
%  ylabel
ylabel({'Axial stress $\sigma_1$ [MPa]'},'Interpreter','latex');
%  xlabel
xlabel({'Axial strain $\varepsilon_1$ [\%]'},'Interpreter','latex');
xlim([-0.006 0])
set(gca,'XDir','reverse','YDir','reverse');





