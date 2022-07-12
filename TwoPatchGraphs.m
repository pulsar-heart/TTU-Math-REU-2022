clear

% birth rates need to be multiplied by initial population

set(0,'DefaultAxesFontSize', 18)
set(gca,'fontsize',18);
tic

% Initial Conditions
i0_w=20;i0_d1=0; i0_h1=0; t0_d1=0; %initial infected
i0_d2=0; i0_h2=0; t0_d2=0;
N0_w=1000; N0_d1=1000; N0_h1=1000; %Initial population
N0_d2=1000; N0_h2=1000;
r0_w=0; r0_d1=0; r0_h1=0; %initial recovered
r0_d2=0; r0_h2=0;
s0_w=N0_w-i0_w; s0_d1=N0_d1-i0_d1; s0_h1=N0_h1-i0_h1; %initial susceptible
s0_d2=N0_d2-i0_d2; s0_h2=N0_h2-i0_h2;


% Parameters
amp_bww=-0.5; beta_ww0 = 0.89;
beta_wd=.51; 
beta_dd=.89; beta_td=.5; beta_th=.207; beta_hh=.078; %transmission rates
gamma_w=52/360; gamma_d=52/360; gamma_h=52/360; %recovery rates
mu_dt=.499; %mutation to human transmissable 
m_w=.00082; m_d=.00082; m_h=.009; %death rates
b_d1=.00082*N0_d1; b_h1=N0_h1*0.0118; %birth rates
b_d2=.00082*N0_d2; b_h2=N0_h2*0.0118;
dt=.001;
%%%
eta = [5 2 1 3]; % Decay rate summer, fall, winter, spring
alpha = 1; % "Environmental infectiousness" % Maybe just here for units?
rho = 10e-3; % Exposure rate
omega = 10e5; % virus shedding rate

% Number of sampe paths, time step, and end time
sim = 10; outbreak=200;
sea=4;
prealloc_size = 52000;
t = [0:dt:dt*prealloc_size];
length(t)
% CTMC Model (Gillespie's algorithm)
for g=1:sea 
   TI=365*g/4;
   totext1(g)=0;
   totext2(g)=0;

for k=1:sim

s_h1 = zeros(1, prealloc_size); s_w = zeros(1,prealloc_size); s_d1 = zeros(1,prealloc_size);
i_h1 = zeros(1,prealloc_size); i_w = zeros(1,prealloc_size); i_d1 = zeros(1,prealloc_size); t_d1 = zeros(1,prealloc_size);
r_h1 = zeros(1,prealloc_size); r_w = zeros(1,prealloc_size); r_d1 = zeros(1,prealloc_size);
s_h2 = zeros(1, prealloc_size);s_d2 = zeros(1,prealloc_size);
i_h2 = zeros(1,prealloc_size); i_d2 = zeros(1,prealloc_size); t_d2 = zeros(1,prealloc_size);
r_h2 = zeros(1,prealloc_size); r_d2 = zeros(1,prealloc_size);
x = zeros(1, prealloc_size);


% Initialize vectors 
t(1)=0; i_w(1)=i0_w; s_w(1)=N0_w-i0_w; r_w(1)=r0_w;
i_d1(1)=i0_d1; i_h1(1)=i0_h1; t_d1(1)=t0_d1; 
s_d1(1)=N0_d1; s_h1(1)=N0_h1;
r_d1(1)=r0_d1; r_h1(1)=r0_h1;
i_d2(1)=i0_d2; i_h2(1)=i0_h2; t_d2(1)=t0_d2; 
s_d2(1)=N0_d2; s_h2(1)=N0_h2;
r_d2(1)=r0_d2; r_h2(1)=r0_h2;
x(1) = 0;
j=1;
N_w = 2; N_d1 = 2; N_h1 = 2;  N_d2 = 2; N_h2 = 2;

  beta_ww= @(j) amp_bww*sin((t(j)+TI)*(3.14/180)) + beta_ww0;
  wd1 = @(j) .5*sin((t(j)+TI)*(3.14/180))+.5;
  wd2 = @(j) -.5*sin((t(j)+TI)*(3.14/180))+.5;

% While loop
while  i_w(j)+i_h1(j)+i_d1(j)+t_d1(j)+i_h2(j)+i_d2(j)+t_d2(j) > 0 && i_w(j)+i_h1(j)+i_d1(j)+t_d1(j)+i_h2(j)+i_d2(j)+t_d2(j) < outbreak && N_w > 1 && N_d1 > 1 && N_h1 > 1 && N_d2 > 1 && N_h2 > 1 
%while j < 5
  u2=rand; % Only need one random variable/ u2=rand; % two uniform random numbers
  


  if sin((t(j)+TI)*(pi/180))-.5>0
       b_w=.0025;
   elseif sin((t(j)+TI)*(pi/180))-.5<=0
       b_w=0;
  end

  N_w=s_w(j)+i_w(j)+r_w(j); % set population sizes
  N_d1=s_d1(j)+i_d1(j)+r_d1(j)+t_d1(j);
  N_h1=s_h1(j)+i_h1(j)+r_h1(j);
  N_d2=s_d2(j)+i_d2(j)+r_d2(j)+t_d2(j);
  N_h2=s_h2(j)+i_h2(j)+r_h2(j);

  ev1=b_w*dt; % birth of wild bird
  ev2=m_w*s_w(j)*dt+ev1;% death of sw
  ev3=m_w*i_w(j)*dt+ev2;% death of Iw
  ev4=dt*(beta_ww(j)*s_w(j)*i_w(j)/N_w + rho*s_w(j)*(1 - exp(- alpha*omega*x(j))))+ev3; % infection of wild bird
  ev5=gamma_w*i_w(j)*dt+ev4;% recovery of wild bird
  ev6=m_w*r_w(j)*dt+ev5;% death of rw
 %Patch 1 Events
  ev7=b_d1*dt+ev6; %birth of chicken 
  ev8=(beta_dd*s_d1(j)*i_d1(j)+beta_wd*wd1(j)*s_d1(j)*i_w(j))/(N_d1)*dt+ev7;% infection of sd by id or iw
  ev9=beta_td*s_d1(j)*t_d1(j)/(N_d1)*dt+ev8;% infection of sd by td
  ev10=gamma_d*i_d1(j)*dt+ev9; % recovery of id
  ev11=gamma_d*t_d1(j)*dt+ev10;% recovery of td
  ev12=m_d*s_d1(j)*dt+ev11;% death of sd
  ev13=m_d*i_d1(j)*dt+ev12; % death of id
  ev14=m_d*t_d1(j)*dt+ev13;% death of td
  ev15=m_d*r_d1(j)*dt+ev14;% death of rd
  ev16=mu_dt*i_d1(j)*dt+ev15; % mutation of id
  ev17=b_h1*dt+ev16;% birth of human
  ev18=(beta_th*s_h1(j)*t_d1(j)+beta_hh*s_h1(j)*i_h1(j))/(N_h1)*dt+ev17;% infection of human
  ev19=m_h*s_h1(j)*dt+ev18; %death of sh
  ev20=m_h*i_h1(j)*dt+ev19;% death of ih
  ev21=m_h*r_h1(j)*dt+ev20;% death of rh
  ev22=gamma_h*i_h1(j)*dt+ev21;
 %Patch 2 Events
  ev23=b_d2*dt+ev22; %birth of chicken 
  ev24=(beta_dd*s_d2(j)*i_d2(j)+beta_wd*wd2(j)*s_d2(j)*i_w(j))/(N_d2)*dt+ev23;% infection of sd by id or iw
  ev25=beta_td*s_d2(j)*t_d2(j)/(N_d2)*dt+ev24;% infection of sd by td
  ev26=gamma_d*i_d2(j)*dt+ev25; % recovery of id
  ev27=gamma_d*t_d2(j)*dt+ev26;% recovery of td
  ev28=m_d*s_d2(j)*dt+ev27;% death of sd
  ev29=m_d*i_d2(j)*dt+ev28; % death of id
  ev30=m_d*t_d2(j)*dt+ev29;% death of td
  ev31=m_d*r_d2(j)*dt+ev30;% death of rd
  ev32=mu_dt*i_d2(j)*dt+ev31; % mutation of id
  ev33=b_h2*dt+ev32;% birth of human
  ev34=(beta_th*s_h2(j)*t_d2(j)+beta_hh*s_h2(j)*i_h2(j))/(N_h2)*dt+ev33;% infection of human
  ev35=m_h*s_h2(j)*dt+ev34; %death of sh
  ev36=m_h*i_h2(j)*dt+ev35;% death of ih
  ev37=m_h*r_h2(j)*dt+ev36;% death of rh
  ev38=gamma_h*i_h2(j)*dt+ev37;

  if u2 <= ev1 % birth of wild bird
    %disp('birth of wild bird')
    s_w(j+1)=s_w(j)+1;
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j);
    i_d1(j+1)=i_d1(j);
    t_d1(j+1)=t_d1(j);
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev1 && u2<=ev2 %death of sw
    %disp('death of sw')
    s_w(j+1)=s_w(j)-1;
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j);
    i_d1(j+1)=i_d1(j);
    t_d1(j+1)=t_d1(j);
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev2 && u2<=ev3 %death of Iw
    %disp('death of Iw')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j)-1;
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j);
    i_d1(j+1)=i_d1(j);
    t_d1(j+1)=t_d1(j);
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev3 && u2<=ev4 %infection of wild bird
    %disp('infection of wild bird')
    s_w(j+1)=s_w(j)-1;
    i_w(j+1)=i_w(j)+1;
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j);
    i_d1(j+1)=i_d1(j);
    t_d1(j+1)=t_d1(j);
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev4 && u2<=ev5 %recovery of wild bird
    %disp('recovery of wild bird')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j)-1;
    r_w(j+1)=r_w(j)+1;
    s_d1(j+1)=s_d1(j);
    i_d1(j+1)=i_d1(j);
    t_d1(j+1)=t_d1(j);
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev5 && u2<=ev6 %death of rw
    %disp('death of rw')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j)-1;
    s_d1(j+1)=s_d1(j);
    i_d1(j+1)=i_d1(j);
    t_d1(j+1)=t_d1(j);
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev6 && u2<=ev7 %birth of chicken
    %disp('birth of chicken')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j)+1;
    i_d1(j+1)=i_d1(j);
    t_d1(j+1)=t_d1(j);
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev7 && u2<=ev8 %infection of sd by id or iw
    %disp('infection of sd by id or iw')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j)-1;
    i_d1(j+1)=i_d1(j)+1;
    t_d1(j+1)=t_d1(j);
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev8 && u2<=ev9 %infection of sd by td
    %disp('infection of sd by td')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j)-1;
    i_d1(j+1)=i_d1(j);
    t_d1(j+1)=t_d1(j)+1;
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev9 && u2<=ev10 %recovery of id
    %disp('recovery of id')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j);
    i_d1(j+1)=i_d1(j)-1;
    t_d1(j+1)=t_d1(j);
    r_d1(j+1)=r_d1(j)+1;
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev10 && u2<=ev11 %recovery of td
    %disp('recovery of td')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j);
    i_d1(j+1)=i_d1(j);
    t_d1(j+1)=t_d1(j)-1;
    r_d1(j+1)=r_d1(j)+1;
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev11 && u2<=ev12 %death of sd
    %disp('death of sd')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j)-1;
    i_d1(j+1)=i_d1(j);
    t_d1(j+1)=t_d1(j);
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev12 && u2<=ev13 %death of id
    %disp('death of id')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j);
    i_d1(j+1)=i_d1(j)-1;
    t_d1(j+1)=t_d1(j);
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev13 && u2<=ev14 %death of td
    %disp('death of td')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j);
    i_d1(j+1)=i_d1(j);
    t_d1(j+1)=t_d1(j)-1;
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev14 && u2<=ev15 %death of rd
    %disp('death of rd')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j);
    i_d1(j+1)=i_d1(j);
    t_d1(j+1)=t_d1(j);
    r_d1(j+1)=r_d1(j)-1;
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev15 && u2<=ev16 %mutation of id
    %disp('mutation of id')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j);
    i_d1(j+1)=i_d1(j)-1;
    t_d1(j+1)=t_d1(j)+1;
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev16 && u2<=ev17 %birth of human
    %disp('birth of human')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j);
    i_d1(j+1)=i_d1(j);
    t_d1(j+1)=t_d1(j);
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j)+1;
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev17 && u2<=ev18 %infection of human
    %disp('infection of human')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j);
    i_d1(j+1)=i_d1(j);
    t_d1(j+1)=t_d1(j);
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j)-1;
    i_h1(j+1)=i_h1(j)+1;
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev18 && u2<=ev19 %death of sh
    %disp('death of sh')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j);
    i_d1(j+1)=i_d1(j);
    t_d1(j+1)=t_d1(j);
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j)-1;
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev19 && u2<=ev20 %death of ih
    %disp('death of ih')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j);
    i_d1(j+1)=i_d1(j);
    t_d1(j+1)=t_d1(j);
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j)-1;
    r_h1(j+1)=r_h1(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev20 && u2<=ev21 %death of rh
    %disp('death of rh')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j);
    i_d1(j+1)=i_d1(j);
    t_d1(j+1)=t_d1(j);
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j);
    r_h1(j+1)=r_h1(j)-1;
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev21 && u2<=ev22 %Recovery of Ih
    %disp('recovery of ih')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d1(j);
    i_d1(j+1)=i_d1(j);
    t_d1(j+1)=t_d1(j);
    r_d1(j+1)=r_d1(j);
    s_h1(j+1)=s_h1(j);
    i_h1(j+1)=i_h1(j)-1;
    r_h1(j+1)=r_h1(j)+1;
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  
  %Patch 2 Events
  elseif u2>ev22 && u2<=ev23 %birth of chicken
    %disp('birth of chicken')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d2(j);
    i_d1(j+1)=i_d2(j);
    t_d1(j+1)=t_d2(j);
    r_d1(j+1)=r_d2(j);
    s_h1(j+1)=s_h2(j);
    i_h1(j+1)=i_h2(j);
    r_h1(j+1)=r_h2(j);
    s_d2(j+1)=s_d2(j)+1;
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev23 && u2<=ev24 %infection of sd by id or iw
    %disp('infection of sd by id or iw')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d2(j);
    i_d1(j+1)=i_d2(j);
    t_d1(j+1)=t_d2(j);
    r_d1(j+1)=r_d2(j);
    s_h1(j+1)=s_h2(j);
    i_h1(j+1)=i_h2(j);
    r_h1(j+1)=r_h2(j);
    s_d2(j+1)=s_d2(j)-1;
    i_d2(j+1)=i_d2(j)+1;
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev24 && u2<=ev25 %infection of sd by td
    %disp('infection of sd by td')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d2(j);
    i_d1(j+1)=i_d2(j);
    t_d1(j+1)=t_d2(j);
    r_d1(j+1)=r_d2(j);
    s_h1(j+1)=s_h2(j);
    i_h1(j+1)=i_h2(j);
    r_h1(j+1)=r_h2(j);
    s_d2(j+1)=s_d2(j)-1;
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j)+1;
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev25 && u2<=ev26 %recovery of id
    %disp('recovery of id')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d2(j);
    i_d1(j+1)=i_d2(j);
    t_d1(j+1)=t_d2(j);
    r_d1(j+1)=r_d2(j);
    s_h1(j+1)=s_h2(j);
    i_h1(j+1)=i_h2(j);
    r_h1(j+1)=r_h2(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j)-1;
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j)+1;
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev26 && u2<=ev27 %recovery of td
    %disp('recovery of td')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d2(j);
    i_d1(j+1)=i_d2(j);
    t_d1(j+1)=t_d2(j);
    r_d1(j+1)=r_d2(j);
    s_h1(j+1)=s_h2(j);
    i_h1(j+1)=i_h2(j);
    r_h1(j+1)=r_h2(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j)-1;
    r_d2(j+1)=r_d2(j)+1;
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev27 && u2<=ev28 %death of sd
    %disp('death of sd')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d2(j);
    i_d1(j+1)=i_d2(j);
    t_d1(j+1)=t_d2(j);
    r_d1(j+1)=r_d2(j);
    s_h1(j+1)=s_h2(j);
    i_h1(j+1)=i_h2(j);
    r_h1(j+1)=r_h2(j);
    s_d2(j+1)=s_d2(j)-1;
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev28 && u2<=ev29 %death of id
    %disp('death of id')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d2(j);
    i_d1(j+1)=i_d2(j);
    t_d1(j+1)=t_d2(j);
    r_d1(j+1)=r_d2(j);
    s_h1(j+1)=s_h2(j);
    i_h1(j+1)=i_h2(j);
    r_h1(j+1)=r_h2(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j)-1;
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev29 && u2<=ev30 %death of td
    %disp('death of td')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d2(j);
    i_d1(j+1)=i_d2(j);
    t_d1(j+1)=t_d2(j);
    r_d1(j+1)=r_d2(j);
    s_h1(j+1)=s_h2(j);
    i_h1(j+1)=i_h2(j);
    r_h1(j+1)=r_h2(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j)-1;
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev30 && u2<=ev31 %death of rd
    %disp('death of rd')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d2(j);
    i_d1(j+1)=i_d2(j);
    t_d1(j+1)=t_d2(j);
    r_d1(j+1)=r_d2(j);
    s_h1(j+1)=s_h2(j);
    i_h1(j+1)=i_h2(j);
    r_h1(j+1)=r_h2(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j)-1;
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev31 && u2<=ev32 %mutation of id
    %disp('mutation of id')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d2(j);
    i_d1(j+1)=i_d2(j);
    t_d1(j+1)=t_d2(j);
    r_d1(j+1)=r_d2(j);
    s_h1(j+1)=s_h2(j);
    i_h1(j+1)=i_h2(j);
    r_h1(j+1)=r_h2(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j)-1;
    t_d2(j+1)=t_d2(j)+1;
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev32 && u2<=ev33 %birth of human
    %disp('birth of human')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d2(j);
    i_d1(j+1)=i_d2(j);
    t_d1(j+1)=t_d2(j);
    r_d1(j+1)=r_d2(j);
    s_h1(j+1)=s_h2(j);
    i_h1(j+1)=i_h2(j);
    r_h1(j+1)=r_h2(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j)+1;
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev33 && u2<=ev34 %infection of human
    %disp('infection of human')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d2(j);
    i_d1(j+1)=i_d2(j);
    t_d1(j+1)=t_d2(j);
    r_d1(j+1)=r_d2(j);
    s_h1(j+1)=s_h2(j);
    i_h1(j+1)=i_h2(j);
    r_h1(j+1)=r_h2(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j)-1;
    i_h2(j+1)=i_h2(j)+1;
    r_h2(j+1)=r_h2(j);
  elseif u2>ev34 && u2<=ev35 %death of sh
    %disp('death of sh')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d2(j);
    i_d1(j+1)=i_d2(j);
    t_d1(j+1)=t_d2(j);
    r_d1(j+1)=r_d2(j);
    s_h1(j+1)=s_h2(j);
    i_h1(j+1)=i_h2(j);
    r_h1(j+1)=r_h2(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j)-1;
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  elseif u2>ev35 && u2<=ev36 %death of ih
    %disp('death of ih')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d2(j);
    i_d1(j+1)=i_d2(j);
    t_d1(j+1)=t_d2(j);
    r_d1(j+1)=r_d2(j);
    s_h1(j+1)=s_h2(j);
    i_h1(j+1)=i_h2(j);
    r_h1(j+1)=r_h2(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j)-1;
    r_h2(j+1)=r_h2(j);
  elseif u2>ev36 && u2<=ev37 %death of rh
    %disp('death of rh')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d2(j);
    i_d1(j+1)=i_d2(j);
    t_d1(j+1)=t_d2(j);
    r_d1(j+1)=r_d2(j);
    s_h1(j+1)=s_h2(j);
    i_h1(j+1)=i_h2(j);
    r_h1(j+1)=r_h2(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j)-1;
  elseif u2>ev37 && u2<=ev38 %Recovery of Ih
    %disp('recovery of ih')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d2(j);
    i_d1(j+1)=i_d2(j);
    t_d1(j+1)=t_d2(j);
    r_d1(j+1)=r_d2(j);
    s_h1(j+1)=s_h2(j);
    i_h1(j+1)=i_h2(j);
    r_h1(j+1)=r_h2(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j)-1;
    r_h2(j+1)=r_h2(j)+1;
  elseif u2>ev38  %No Event
    %disp('recovery of ih')
    s_w(j+1)=s_w(j);
    i_w(j+1)=i_w(j);
    r_w(j+1)=r_w(j);
    s_d1(j+1)=s_d2(j);
    i_d1(j+1)=i_d2(j);
    t_d1(j+1)=t_d2(j);
    r_d1(j+1)=r_d2(j);
    s_h1(j+1)=s_h2(j);
    i_h1(j+1)=i_h2(j);
    r_h1(j+1)=r_h2(j);
    s_d2(j+1)=s_d2(j);
    i_d2(j+1)=i_d2(j);
    t_d2(j+1)=t_d2(j);
    r_d2(j+1)=r_d2(j);
    s_h2(j+1)=s_h2(j);
    i_h2(j+1)=i_h2(j);
    r_h2(j+1)=r_h2(j);
  end
  x(j + 1) = x(j)*exp(-eta(g)*dt) + i_w(1, j)*(1 - exp(-eta(g)*dt))/eta(g);
  j=j+1;
  
end  

% resize vectors to remove zeros at end
s_h1 = s_h1( 1:j); s_w = s_w( 1:j); s_d1 = s_d1( 1:j); 
i_h1 = i_h1( 1:j); i_w = i_w( 1:j); i_d1 = i_d1( 1:j); t_d1 = t_d1( 1:j); 
r_h1 = r_h1( 1:j); r_w = r_w( 1:j); r_d1 = r_d1( 1:j); 

s_h2 = s_h2( 1:j); s_d2 = s_d2( 1:j); 
i_h2 = i_h2( 1:j); i_d2 = i_d2( 1:j); t_d2 = t_d2( 1:j); 
r_h2 = r_h2( 1:j); r_d2 = r_d2( 1:j); 
t_graph = t(1:j);


figure(g)
if k<=7
    stairs(t_graph,i_d2,'color',rand(1,3),'Linewidth',2);
    hold on

end

if i_h1(j)+i_w(j)+i_d1(j)+t_d1(j)==0
    totext1(g)=totext1(g)+1;
elseif i_h2(j)+i_w(j)+i_d2(j)+t_d2(j)==0
    totext2(g)=totext2(g)+1;
end
end
end
figure(1)
title('CTMC Summer with Periodicity')
ylabel('I_w(t)');
xlabel('Time (days)')
hold off

figure(2)
title('CTMC Fall with Periodicity')
ylabel('I_w(t)');
xlabel('Time (days)')
hold off

figure(3)
title('CTMC Winter with Periodicity')
ylabel('I_w(t)');
xlabel('Time (days)')
hold off

figure(4)
title('CTMC Spring with Periodicity')
ylabel('I_w(t)');
xlabel('Time (days)')
hold off

sim_prob_ext_summer_1=totext1(1)/sim
sim_prob_ext_fall_1=totext1(2)/sim
sim_prob_ext_winter_1=totext1(3)/sim
sim_prob_ext_spring_1=totext1(4)/sim

sim_prob_ext_summer_2=totext2(1)/sim
sim_prob_ext_fall_2=totext2(2)/sim
sim_prob_ext_winter_2=totext2(3)/sim
sim_prob_ext_spring_2=totext2(4)/sim

R0h=beta_hh/(gamma_h+m_h)
R0w=beta_ww0/(gamma_w+m_w)
R0t=beta_dd/(gamma_d+m_d)
R0d=beta_dd/(gamma_d+m_d+mu_dt)

toc

