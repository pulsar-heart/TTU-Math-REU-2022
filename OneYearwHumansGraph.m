%%%%%%%%%%%%%%%%%%%%%%%%%% 
% CTMC Model Environmental Transmission Full Stochastic
% %%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%dbclear all

% time (x*360 --> x=number of years)
  dt=10e-4; tmax=360; nyears = 1;

% birth rates need to be multiplied by initial population
% Initial conditions
N0_w1=2500; N0_w2=2500; N0_d=2000; N0_h=1000; %Initial population
N_w1 = N0_w1; N_w2 = N0_w2; N_d = N0_d; N_h = N0_h;

i0_w1=1;i0_w2=0;i0_d=0; i0_h=0; t0_d=0; %initial infected
r0_w1=0; r0_w2=0; r0_d=0; r0_h=0; %initial recovered
s0_w1=N0_w1-i0_w2; s0_w2=N0_w2-i0_w2; 
s0_d=N0_d-i0_d-t0_d; s0_h=N0_h-i0_h; %initial susceptible

% Parameters
% R0w=0.5; R0d=0.5; R0t=0.5; R0h=0.5;
% recovery rates
  gamma_w=1/11; gamma_d=1/11; gamma_h=1/11; 
% death rates
  m_w=1/(3*360); m_d=1/360; m_h=1/(76*360);
% birth rates
  b_d=50/9; b_h = m_h*N0_h; b_w0 = 4.58703;
% transmission
  beta_ww0 = 0.046;  beta_dd=0.046; beta_td=0.046;  beta_hh = 3*0.046;
  %beta_hh=4.6*10^(-5);
  beta_wd=4.6*10^(-3);  beta_th=0.0046; %4.6*10^(-4); %contact based
% environmental transmission
  rho_w = 10^(-3)/360; rho_wd = 5*10^(-4)/360; 
  rho_d = 10^(-3)/360;  rho_h = 10^(-3)/360;
  omega_w = 1; omega_d = 1; omega_t = 1; 
% mutation to human transmissible
  mu_dt=2*0.0005; 

  first_extinct = [0 0 0 0];

  % count when j==0;
% Number of sampe paths, time step, and end time
sim = 1; outbreak=500;
sea=4;
t = [0:dt:tmax];

% CTMC Model (Gillespie's algorithm)

for g=1:sea
%  g = 1;
%  TI = 0;

TI=90*(g - 1) + 60;

for k=1:10

    % prealloc
complete_season = 0;
prealloc_size = tmax/dt;
s_w1=zeros(1, prealloc_size); i_w1=zeros(1, prealloc_size); r_w1=zeros(1, prealloc_size);
s_w2=zeros(1, prealloc_size); i_w2=zeros(1, prealloc_size); r_w2=zeros(1, prealloc_size); 
s_d=zeros(1, prealloc_size); i_d=zeros(1, prealloc_size); t_d=zeros(1, prealloc_size); r_d=zeros(1, prealloc_size);
s_h=zeros(1, prealloc_size); i_h=zeros(1, prealloc_size); r_h=zeros(1, prealloc_size);
v_n=zeros(1, prealloc_size); v_s=zeros(1, prealloc_size); v_d=zeros(1, prealloc_size); v_t=zeros(1, prealloc_size);
% Initialize vectors
beta_ww(1)=beta_ww0; 
s_w1(1)=s0_w1; i_w1(1)=i0_w1; r_w1(1)=r0_w1; N_w1=s_w1(1)+i_w1(1)+r_w1(1);
s_w2(1)=s0_w2; i_w2(1)=i0_w2; r_w2(1)=r0_w2; N_w2=s_w2(1)+i_w2(1)+r_w2(1);
s_d(1)=s0_d; i_d(1)=i0_d; t_d(1)=t0_d; r_d(1)=r0_d; N_d=s_d(1)+i_d(1)+t_d(1)+r_d(1);
s_h(1)=s0_h; i_h(1)=i0_h; r_h(1)=r0_h; n_h=s_h(1)+i_h(1)+r_h(1);
 v_n(1)=0;  v_s(1)=0; v_d(1)=0; v_t(1)=0;
j = 1;
infected = 1;
% While loop
%while  infected > 0 && infected  < outbreak && N_w1 > 1 && N_w2 > 1&& N_d > 1 && N_h > 1 && j < tmax/dt
while 1
  u2=rand; % Only need one random variable/ u2=rand; % two uniform random numbers
   time = t(j) + TI;
   b_w = b_w0*(tanh(57*sin(pi/180*(time - 106)) - 40) + 1);
   beta_ww = 0.046*(1 + 0.25*sin(pi/180*(time - 180)));
   m = 0.5*(tanh(180*sin(pi/180*(time - 234)) - 4) + 1);
   eta_s = 4/35*(1 + 3/4*sin(pi/180*(time - 90)));
   eta_n = 11/525*(1 + 4/11*sin(pi/180*(time - 90)));

  N_w1=s_w1(j) +i_w1(j) +r_w1(j) ; % set population sizes
  N_w2=s_w2(j) +i_w2(j) +r_w2(j) ; 
  N_d=s_d(j) +i_d(j) +r_d(j) +t_d(j) ;
  N_h= s_h(j) + i_h(j) +r_h(j) ;

  ev1=b_w*dt; % birth of wild bird 1
  ev2=b_w*dt+ev1; % birth wild patch 2
  ev3=m_w*s_w1(j) *dt+ev2;% death of sw 1
  ev4=m_w*s_w2(j) *dt+ev3;% death of sw 2
  ev5=m_w*i_w1(j) *dt+ev4;% death of Iw 1
  ev6=m_w*i_w2(j) *dt+ev5;% death of Iw 2
  ev7=dt*(beta_ww *s_w1(j) *i_w1(j) /N_w1 + m*beta_ww *s_w1(j) *i_w2(j) /N_w1  + ...
      (1-m) *rho_w*s_w1(j) *v_n(j)  + m*rho_w*s_w1(j) *v_s(j) )+ev6; % infection of wild bird 1
  ev8=dt*(beta_ww *s_w2(j) *i_w2(j) /N_w2 + m*beta_ww *s_w2(j) *i_w1(j) /N_w2  + ...
      rho_w*s_w2(j) *v_s(j) )+ev7; % infection of wild bird  2
  ev9=gamma_w*i_w1(j)*dt+ev8; % recovery of iw1
  ev10=gamma_w*i_w2(j)*dt+ev9; % recovery of iw2
  ev11=m_w*r_w1(j) *dt+ev10;% death of rw1
  ev12=m_w*r_w2(j) *dt+ev11;% death of rw2
  ev13=b_d*dt+ev12; %birth of chicken
  ev14=((beta_dd*s_d(j) *i_d(j) +m*beta_wd*s_d(j) *i_w1(j) +beta_wd*s_d(j) *i_w2(j) )/N_d + ...
     s_d(j) *(rho_d*v_d(j) +rho_wd*v_s(j) ))*dt+ev13;% infection of sd by id or iw % Need Rho_wd
  ev15=(beta_td*s_d(j) *t_d(j) /N_d + rho_d*s_d(j) *v_t(j) )*dt+ev14;% infection of sd by td
  ev16=gamma_d*i_d(j) *dt+ev15; % recovery of id
  ev17=gamma_d*t_d(j) *dt+ev16;% recovery of td
  ev18=m_d*s_d(j) *dt+ev17;% death of sd
  ev19=m_d*i_d(j) *dt+ev18; % death of id
  ev20=m_d*t_d(j) *dt+ev19;% death of td
  ev21=m_d*r_d(j) *dt+ev20;% death of rd
  ev22=mu_dt*i_d(j) *dt+ev21; % mutation of id
  ev23=b_h*dt+ev22;% birth of human
  ev24=((beta_th*s_h(j) *t_d(j) +beta_hh*s_h(j) *i_h(j) )/N_h + rho_h*s_h(j) *v_t(j) )*dt+ev23;% infection of human
  ev25=m_h*s_h(j) *dt+ev24; %death of sh
  ev26=m_h*i_h(j) *dt+ev25;% death of ih
  ev27=m_h*r_h(j) *dt+ev26;% death of rh
  ev28=gamma_h*i_h(j) *dt+ev27; % recovery of human
  ev29=(1-m)*omega_w*i_w1(j) *dt+ev28; % shed infectious dose v_w 1
  ev30=(m*omega_w*i_w1(j)  + omega_w*i_w2(j) )*dt+ev29; % shed infectious dose v_w 2
  ev31=omega_d*i_d(j) *dt+ev30; % shed infectious dose v_d
  ev32=omega_t*t_d(j) *dt+ev31; % shed infectious dose v_t
  ev33=eta_n*v_n(j) *dt+ev32; % decay infectious dose v_n
  ev34=eta_s*v_s(j) *dt+ev33; % decay infectious dose v_s
  ev35=eta_s*v_d(j) *dt+ev34; % decay infectious dose v_d
  ev36=eta_s*v_t(j) *dt+ev35; % decay infectious dose v_t

  if u2 <= ev1 % birth of wild bird 1
    %disp('birth of wild bird  1')
    s_w1(j+1) =s_w1(j) +1;
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);



  elseif u2 <= ev2 % birth of wild bird2
    %disp('birth of wild bird2')
    s_w2(j+1) =s_w2(j) +1;
    s_w1(j+1) = s_w1(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
     v_n(j+1)= v_n(j);
     v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev3  %death of sw1
    %disp('death of sw1')
    s_w1(j+1) =s_w1(j) -1;
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
     v_n(j+1)= v_n(j);
     v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev4  %death of sw2
    %disp('death of sw2')
    s_w2(j+1) =s_w2(j) -1;
    s_w1(j+1) = s_w1(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
     v_n(j+1)= v_n(j);
     v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev5 %death of Iw1
    %disp('death of Iw1')
    i_w1(j+1) =i_w1(j) -1;
    s_w1(j+1) = s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
     v_n(j+1)= v_n(j);
     v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev6 %death of Iw2
    %disp('death of Iw2')
    i_w2(j+1) =i_w2(j) -1;
    s_w1(j+1) = s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);

    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
     v_n(j+1)= v_n(j);
     v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);
    

  elseif u2<=ev7  %infection of wild bird1
    %disp('infection of wild bird')
    s_w1(j+1) =s_w1(j) -1;
    i_w1(j+1) =i_w1(j) +1;

    s_w2(j+1) = s_w2(j);

    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
     v_n(j+1)= v_n(j);
     v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev8  %infection of wild bird2
    %disp('infection of wild bird')
    s_w2(j+1) =s_w2(j) -1;
    i_w2(j+1) =i_w2(j) +1;
    s_w1(j+1) = s_w1(j);

    i_w1(j+1)=i_w1(j);

    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
     v_n(j+1)= v_n(j);
     v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev9  %recovery of wild bird
    %disp('recovery of wild bird')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j) - 1;
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j) + 1;
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev10  %recovery of wild bird
    %disp('recovery of wild bird')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j) - 1;
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j) + 1;
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);
  elseif u2<=ev11 %death of rw
    %disp('death of rw')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j) - 1;
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev12  %death of rw
    %disp('death of rw')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j) - 1;
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev13  %birth of chicken
    %disp('birth of chicken')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j) + 1;
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev14  %infection of sd by id or iw
    %disp('infection of sd by id or iw')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j) - 1;
    i_d(j+1)=i_d(j) + 1;
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev15  %infection of sd by td
    %disp('infection of sd by td')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j) - 1;
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j) + 1;
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev16  %recovery of id
    %disp('recovery of id')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j) - 1;
    r_d(j+1)=r_d(j) + 1;
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev17  %recovery of td
    %disp('recovery of td')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j) + 1;
    t_d(j+1)=t_d(j) - 1;
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev18 %death of sd
    %disp('death of sd')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j) - 1;
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev19  %death of id
    %disp('death of id')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j) - 1;
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev20  %death of td
    %disp('death of td')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j) - 1;
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev21  %death of rd
    %disp('death of rd')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j) - 1;
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev22  %mutation of id
    %disp('mutation of id')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j) - 1;
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j) + 1;
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev23  %birth of human
    %disp('birth of human')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j) + 1;
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev24  %infection of human
    %disp('infection of human')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j) - 1;
    i_h(j+1)=i_h(j) + 1;
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev25   %death of sh
    %disp('death of sh')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j) - 1;
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev26 %death of ih
    %disp('death of ih')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j) - 1;
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev27  %death of rh
    %disp('death of rh')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j) - 1;
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev28  %Recovery of Ih
    %disp('recovery of ih')
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j) - 1;
    r_h(j+1)=r_h(j) + 1;
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev29 % shed infectious dose  v_n
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j) + 1;
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev30 % shed infectious dose  v_s
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j) + 1;
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev31 % shed infectious dose v_d
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j) + 1;
    v_t(j+1)=v_t(j);

  elseif u2<=ev32 % shed infectious dose v_t
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j) + 1;

  elseif u2<=ev33 % decay infectious dose  v_n
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j) - 1;
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev34 % decay infectious dose  v_s
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j) - 1;
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);

  elseif u2<=ev35 % decay infectious dose v_d
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j) - 1;
    v_t(j+1)=v_t(j);

  elseif u2<=ev36 % decay infectious dose v_t
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j) - 1;
  else
    s_w1(j+1) =s_w1(j);
    s_w2(j+1) = s_w2(j);
    i_w1(j+1)=i_w1(j);
    i_w2(j+1)=i_w2(j);
    r_w1(j+1)=r_w1(j);
    r_w2(j+1)=r_w2(j);
    s_d(j+1)=s_d(j);
    i_d(j+1)=i_d(j);
    r_d(j+1)=r_d(j);
    t_d(j+1)=t_d(j);
    s_h(j+1)=s_h(j);
    i_h(j+1)=i_h(j);
    r_h(j+1)=r_h(j);
    v_n(j+1)= v_n(j);
    v_s(j+1)= v_s(j);
    v_d(j+1)=v_d(j);
    v_t(j+1)=v_t(j);



  end
  infected = i_w1(j) + i_w2(j) + i_d(j) + t_d(j) + i_h(j) +  v_n(j) +  v_s(j) + v_d(j) + v_t(j);
  if infected > outbreak
      complete_season = 1;
  end
     
  if j >= tmax/dt
     break
  end %end if block   

  j = j+1;
end  % end while
if complete_season == 1
 
    % Trim vectors
    s_w1 =s_w1(1:j);
    s_w2 = s_w2(1:j);
    i_w1=i_w1(1:j);
    i_w2=i_w2(1:j);
    r_w1=r_w1(1:j);
    r_w2=r_w2(1:j);
    s_d=s_d(1:j);
    i_d=i_d(1:j);
    r_d=r_d(1:j);
    t_d=t_d(1:j);
    s_h=s_h(1:j);
    i_h=i_h(1:j);
    r_h=r_h(1:j);
    v_n= v_n(1:j);
    v_s= v_s(1:j);
    v_d=v_d(1:j);
    v_t=v_t(1:j);
    yrs = t(1:j)/360;
    v_scales = v_s*2*10^(-2);
    v_scalen = v_n*2*10^(-2);
    subplot(2, 2, g);
    hold on
    TIyrs=(90*(g - 1) + 60)/360;
    yyaxis left
    set(gca, 'YColor', 'k')
    stairs(yrs + TIyrs, i_d, 'Color', 'black', 'Linestyle', '--', 'Linewidth',2);
    stairs(yrs + TIyrs, i_h, 'Color', 'black', 'Linestyle', '-', 'Linewidth', 2);
    yyaxis right
    ylim([0 15])
    xlim([TIyrs, 1 + TIyrs]) 
    stairs(yrs + TIyrs, t_d, 'Color', 'red', 'Linewidth',2);
    complete_season = 0;
    break

end % end if for plotting

end % end for loop of sims
end % end seasonal loop

figure(1)
subplot(2, 2, 1)
title('CTMC Spring with Periodicity')
legend('I_d', 'I_h', 'T_d')
xlabel('years')
set(gca, 'box', 'on')

subplot(2, 2, 2)
title('CTMC Summer with Periodicity')
legend('I_d', 'I_h', 'T_d')
xlabel('years')
set(gca, 'box', 'on')

subplot(2, 2, 3)
title('CTMC Fall with Periodicity')
legend('I_d', 'I_h', 'T_d')
xlabel('years')
set(gca, 'box', 'on')

subplot(2, 2, 4)
title('CTMC Winter with Periodicity')
legend('I_d', 'I_h', 'T_d')
xlabel('years')
set(gca, 'box', 'on')
hold off


