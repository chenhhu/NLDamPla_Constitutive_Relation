% This code is developed to give the analytical solution of the 
% nonlocal continuum damage-plasticity model, see Chen et al., 2024 (Journal). 
% Stress state: Uniaxial compression stree state
% test data-Yang et al. 2008
clear all
close all
clc
%%
% csvread(filename,R1,C1,[R1 C1 R2 C2]) extracted from experimental datas 
y1=csvread('pc0MPa.csv',0,0,[0,0,49,1])';
y2=csvread('pc5MPa.csv',0,0,[0,0,42,1])';
y3=csvread('pc10MPa.csv',0,0,[0,0,41,1])'; 
y4=csvread('pc20MPa.csv',0,0,[0,0,40,1])';
%
color=parula(4);
plot(y1(1,:),y1(2,:),'o','color',color(1,:));
hold on
plot(y2(1,:),y2(2,:),'o','color',color(2,:));
plot(y3(1,:),y3(2,:),'o','color',color(3,:));
plot(y4(1,:),y4(2,:),'o','color',color(4,:));

%%
syms chi_p E c_f alpha_p c_0 h_p beta kappa kappa_i sigma1 A B D sigma2
% local equivalent strain
kappa_fun0=sqrt( 2*chi_p/E*( c_f*alpha_p+(c_f-c_0)/h_p*(exp(-h_p*alpha_p)-1) ) );
kappa_fun1=str2func(['@(chi_p,E,c_f,alpha_p,c_0,h_p)',vectorize(kappa_fun0)]); % Construct function handle from function 
% damage evolution
D_fun0=1-exp(-beta*(kappa-kappa_i));
D_fun1=str2func(['@(beta,kappa,kappa_i)',vectorize(D_fun0)]); % Construct function handle from function 
% The following is the yield function under Uniaxial compression
f_fun0 = sqrt((sigma1-sigma2)^2/3)+A/3*(sigma1+2*sigma2)-B*(1-chi_p*D)*(c_f-(c_f-c_0)*exp(-h_p*alpha_p))==0; 
sigma1_fun0=solve(f_fun0,sigma1); % sigma1 = f(alpha_p)
sigma1_fun1=str2func(['@(A,B,chi_p,D,c_f,c_0,h_p,alpha_p,sigma2)',vectorize(sigma1_fun0(1))]); % Construct function handle from function "d2sigma1".
%% Initiation of material parameters. 
E=66000; %[MPa]
v=0.2;
c_0=16.6875;
theta=pi/180*40;
A=6*sin(theta)/sqrt(3)/(3+sin(theta));
B=6*cos(theta)/sqrt(3)/(3+sin(theta));
C=6*sin(theta)/sqrt(3)/(3+sin(theta));

sigma_c=3;                        
kappa_i=sigma_c/E;
c_f=193.5139;
h_p=191.0201;
beta=494.8072;

sigma2_list=[0 -5 -10 -20];
chi_p_list=[0.995 0.9167 0.8641 0.8382];

for j=1:1:4
sigma2=sigma2_list(j);
chi_p=chi_p_list(j);
%% Load by increasing d
i=0;
for alpha_p=0.000001:0.0002:0.035
    i=i+1;
    kappa=kappa_fun1(chi_p,E,c_f,alpha_p,c_0,h_p);
    if kappa<kappa_i
        D(i)=0;
    else
        D(i)=D_fun1(beta,kappa,kappa_i);
    end
    sigma1(i)=sigma1_fun1(A,B,chi_p,D(i),c_f,c_0,h_p,alpha_p,sigma2); % Axial stress
    N=-sqrt(3)/3+C/3; % Axial plastic flow direction
    strain1(i)=alpha_p/B*N+sigma1(i)/E-2*v*sigma2/E; % Axial strain   
end  
strain1_0=sigma2/E-2*v*sigma2/E;
plot([0,strain1-strain1_0],[0,sigma1-sigma2],'-','color',color(j,:));
end

%% 
%  ylabel
ylabel({'Axial stress $\sigma_1$ [MPa]'},'Interpreter','latex');
%  xlabel
xlabel({'Axial strain $\varepsilon_1$ [\%]'},'Interpreter','latex');
xlim([-0.02 0])
set(gca,'XDir','reverse','YDir','reverse');



