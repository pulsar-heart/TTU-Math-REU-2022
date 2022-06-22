%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Royce CTMC Model
% %%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

% birth rates need to be multiplied by initial population

set(0,'DefaultAxesFontSize', 18)
set(gca,'fontsize',18);

% Parameters
beta_ww=.89; beta_wd=.51; beta_dd=.89; beta_td=.5; beta_th=.5; beta_hh=.078; %transmission rates
gamma_w=.981; gamma_d=.981; gamma_h=0.091; %recovery rates
mu_dt=.499; %mutation to human transmissable 
m_w=1; m_d=1; m_h=.009; %death rates
b_w=1000; b_d=1000; b_h=1000*0.0118; %birth rates

% Initial Conditions
i0_w=2;i0_d=0; i0_h=0; t0_d=0; %initial infected
N0_w=1000;N0_d=1000; N0_h=1000; %Initial population
r0_w=0; r0_d=0; r0_h=0; %initial recovered
s0_w=1000; s0_d=N0_d; s0_h=N0_h; %initial susceptible

totext=0; % count when j==0;
% Number of sampe paths, time step, and end time
sim = 10; outbreak=50;



% CTMC Model (Gillespie's algorithm)
for k=1:sim
clear s_h s_w s_d i_w i_d i_h t_d r_h r_w r_d t N_h N_d N_w j
prealloc_size = 5000;
s_h = zeros(prealloc_size, 1); s_w = zeros(prealloc_size,1); s_d = zeros(prealloc_size,1);
i_h = zeros(prealloc_size, 1); i_w = zeros(prealloc_size,1); i_d = zeros(prealloc_size,1); t_d = zeros(prealloc_size,1);
r_h = zeros(prealloc_size, 1); r_w = zeros(prealloc_size,1); r_d = zeros(prealloc_size,1);
t = zeros(prealloc_size, 1);

% Initialize vectors 
t(1)=0; i_w(1)=i0_w; i_d(1)=i0_d; i_h(1)=i0_h; t_d(1)=t0_d; 
s_w(1)=N0_w-i0_w; s_d(1)=N0_d; s_h(1)=N0_h;
r_w(1)=r0_w; r_d(1)=r0_d; r_h(1)=r0_h;
j=1;
N_w = 2; N_d = 2; N_h = 2; 

% While loop
while  i_h(j)+i_d(j)+i_w(j)+t_d(j) > 0 && i_h(j)+i_d(j)+i_w(j)+t_d(j) < outbreak && N_w > 1 && N_d > 1 && N_h > 1 

  u1=rand;  u2=rand; % two uniform random numbers

  N_w=s_w(j)+i_w(j)+r_w(j); % set population sizes
  N_d=s_d(j)+i_d(j)+r_d(j)+t_d(j);
  N_h=s_h(j)+i_h(j)+r_h(j);

  sum = b_w + m_w*s_w(j) + m_w*i_w(j) + beta_ww*s_w(j)*i_w(j)/N_w + gamma_w*i_w(j) + m_w*r_w(j) + ...
      b_d + (beta_dd*s_d(j)*i_d(j)+beta_wd*s_d(j)*i_w(j))/N_d + beta_td*s_d(j)*t_d(j)/N_d + gamma_d*i_d(j) + ...
      gamma_d*t_d(j) + m_d*s_d(j) + m_d*i_d(j) + m_d*t_d(j) + m_d*r_d(j) + mu_dt*i_d(j) + b_h + ...
      (beta_th*s_h(j)*t_d(j)+beta_hh*s_h(j)*i_h(j))/N_h + m_h*s_h(j) + m_h*i_h(j) + m_h*r_h(j) + gamma_h*i_h(j);
  ev1=b_w/sum; % birth of wild bird
  ev2=m_w*s_w(j)/sum+ev1;% death of sw
  ev3=m_w*i_w(j)/sum+ev2;% death of Iw
  ev4=beta_ww*s_w(j)*i_w(j)/(N_w*sum)+ev3; % infection of wild bird
  ev5=gamma_w*i_w(j)/sum+ev4;% recovery of wild bird
  ev6=m_w*r_w(j)/sum+ev5;% death of rw
  ev7=b_d/sum+ev6; %birth of chicken
  ev8=(beta_dd*s_d(j)*i_d(j)+beta_wd*s_d(j)*i_w(j))/(N_d*sum)+ev7;% infection of sd by id or iw
  ev9=beta_td*s_d(j)*t_d(j)/(N_d*sum)+ev8;% infection of sd by td
  ev10=gamma_d*i_d(j)/sum+ev9; % recovery of id
  ev11=gamma_d*t_d(j)/sum+ev10;% recovery of td
  ev12=m_d*s_d(j)/sum+ev11;% death of sd
  ev13=m_d*i_d(j)/sum+ev12; % death of id
  ev14=m_d*t_d(j)/sum+ev13;% death of td
  ev15=m_d*r_d(j)/sum+ev14;% death of rd
  ev16=mu_dt*i_d(j)/sum+ev15; % mutation of id
  ev17=b_h/sum+ev16;% birth of human
  ev18=(beta_th*s_h(j)*t_d(j)+beta_hh*s_h(j)*i_h(j))/(N_h*sum)+ev17;% infection of human
  ev19=m_h*s_h(j)/sum+ev18; %death of sh
  ev20=m_h*i_h(j)/sum+ev19;% death of ih
  ev21=m_h*r_h(j)/sum+ev20;% death of rh
  ev22=gamma_h*i_h(j)/sum+ev21;
  
  

  t(j+1)=t(j)-log(u1)/sum; % interevent time

  if u2 <= ev1 % birth of wild bird
    s_w(j+1)=s_w(j)+1;
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev1 && u2<=ev2 %death of sw
    s_w(j+1)=s_w(j)-1;
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev2 && u2<=ev3 %death of Iw
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j)-1;
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev3 && u2<=ev4 %infection of wild bird
    s_w(j+1)=s_w(j)-1;
    i_w(j+1)=i_w(j)+1;
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev4 && u2<=ev5 %recovery of wild bird
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j)-1;
    r_w(j+1)=r_w(j)+1;
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev5 && u2<=ev6 %death of rw
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j)-1;
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev6 && u2<=ev7 %birth of chicken
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j)+1;
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev7 && u2<=ev8 %infection of sd by id or iw
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j)+1;
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev8 && u2<=ev9 %infection of sd by td
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j)+1;
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    elseif u2>ev9 && u2<=ev10 %recovery of id
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j)-1;
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j)+1;
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev10 && u2<=ev11 %recovery of td
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j)-1;
    r_d(j+1)=r_d(j)+1;
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev11 && u2<=ev12 %death of sd
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j)-1;
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    elseif u2>ev12 && u2<=ev13 %death of id
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j)-1;
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev13 && u2<=ev14 %death of td
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j)-1;
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev14 && u2<=ev15 %death of rd
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j)-1;
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
     elseif u2>ev15 && u2<=ev16 %mutation of id
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j)-1;
    t_d(j+1)=t_d(j)+1;
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev16 && u2<=ev17 %birth of human
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j)+1;
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev17 && u2<=ev18 %infection of human
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j)+1;
    r_h(j+1)=r_h(j);
     elseif u2>ev18 && u2<=ev19 %death of sh
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j)-1;
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
  elseif u2>ev19 && u2<=ev20 %death of ih
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j)-1;
    r_h(j+1)=r_h(j);
  elseif u2>ev20 && u2<=ev21 %death of rh
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j)-1;
  elseif u2>ev21 && u2<=ev22 %Recovery of Ih
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    t_d(j+1)=t_d(j);
    r_d(j+1)=r_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j)-1;
    r_h(j+1)=r_h(j)+1;
  end
  j=j+1;
  
end  
% resize vectors to remove zeros at end
s_h = s_h(1:j); s_w = (1:j); s_d = (1:j); 
i_h = (1:j); i_w = (1:j); i_d = (1:j); t_d = (1:j); 
r_h = (1:j); r_w = (1:j); r_d = (1:j); 
t = (1:j); 

if k==1
    disp(size(t))
stairs(t,i_h,'color',rand(1,3),'Linewidth',2);
hold on
end
if i_h(j)+i_w(j)+i_d(j)==0
    totext=totext+1;
end
end

plot(t,i_h,'k--','LineWidth',2);
title('CTMC Royce Model')
ylabel('I_h(t)');
xlabel('Time (days)')
hold off

