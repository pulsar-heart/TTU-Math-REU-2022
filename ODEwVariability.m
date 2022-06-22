%%%%%%%%%%%%%%%%%%%%%%%%%% 
% SIR Epidemic Model
% with disease-relatd deaths
% Sample Paths for CTMC  models
% ODE model:
% ds/dt=-beta is/N;
% di/dt= beta is/N-(gam +alp)*i;
% dr/dt= gam*i;
% %%%%%%%%%%%%%%%%%%%%%%%%
dbstop if error
clear all
close all
set(0,'DefaultAxesFontSize', 18)
set(gca,'fontsize',18);
%Parameters and Initial Condition

beta_ww0 = 2;
beta_wd=.51; beta_dd=.89; beta_td=.89; beta_th=.207; beta_hh=0.078; %transmission rates
gamma_w=.2; gamma_d=.981; gamma_h=.091; %recovery rates
mu_dt=.499; %mutation to human transmissable 
m_w=.0041; m_d=1; m_h=.009; %death rates
 b_d=1000; b_h=.0118; %birth rates
i0_w=500;i0_d=0; i0_h=0; t0_d=0; %initial infected
N0_w=100000; N0_d=100000; N0_h=100000; %Initial population
r0_w=0; r0_d=0; r0_h=0; %initial recovered
s0_w=500; s0_d=1000; s0_h=1000; %initial susceptible


totext=0; % count when i==0;
%Number of sampe paths, time step, and end time
sim=100; 
dt=0.01; time=800;  outbreak=50;
%sim=10^5; %use for calc of dis extinction
% ODE Model (Euler's Method)

beta_ww_func = @(k) 0.25*sin(k*dt*pi/182) + beta_ww0;

q(1)=s0_w; s(1)=i0_w; a(1)=r0_w; Nw(1)=q(1)+s(1)+a(1);
plot_vec(1) = 0;
for k=1:time/dt
   if sin(k*dt*(pi/180))-.5>0
       b_w=.0123;
       plot_vec(k+1) = 500;
   elseif sin(k*dt*(pi/180))-.5<=0
       b_w=0;
       plot_vec(k+1) = 0;
   end
  w1=Nw(k)*b_w-beta_ww_func(k)*q(k)*s(k)/Nw(k)-m_w*q(k);
  w2=beta_ww_func(k)*q(k)*s(k)/Nw(k)-(gamma_w+m_w)*s(k);
  w3=gamma_w*s(k)-m_w*a(k);
  q(k+1)=q(k)+dt*w1;
  s(k+1)=s(k)+dt*w2;
  a(k+1)=a(k)+dt*w3;
  Nw(k+1)=q(k+1)+s(k+1)+a(k+1);
end

plot([0:dt:time],q,'k--',[0:dt:time],s,'b-.',[0:dt:time],a, 'g--', [0:dt:time], plot_vec,'r--','LineWidth',2);
axis([0,time,0,2000]);
title('ODE SIR with periodic birth rate')
legend('susceptible','infected','recovered')
ylabel('wild birds');
xlabel('Time (days)')
hold off