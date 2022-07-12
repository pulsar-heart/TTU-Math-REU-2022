%%%%%%%%%%%%%%%%%%%%%%%%%% 
% CTMC Model Environmental Transmission Full Stochastic FAST
% %%%%%%%%%%%%%%%%%%%%%%%%
clear
dbstop if error
tic
% birth rates need to be multiplied by initial population

set(0,'DefaultAxesFontSize', 18);
set(gca,'fontsize',18);

N0_w=1000;N0_d=1000; N0_h=1000; %Initial population


R0w=2.5; R0d=2.5; R0t=5; R0h=2.5;
gamma_w=1/11; gamma_d=1/11; gamma_h=1/11; %recovery rates
mu_dt=.499; %mutation to human transmissible 
m_w=1/(3*360); m_d=1/360; m_h=1/(70*360); %death rates
beta_ww0 = R0w*(m_w+gamma_w); amp_bww=0.25*beta_ww0; %seasonal transmission
beta_dd=R0d*(m_w+gamma_w+mu_dt); beta_td=R0t*(m_d+gamma_d); beta_hh=R0h*(m_h+gamma_h); %R0 transmission rates
beta_wd=0.25*beta_dd;  beta_th=0.5*beta_td; %contact based transmission
b_w=1/(3*360)*N0_w; b_d=m_d*N0_d; b_h=N0_h*m_h; %birth rates
dt=10e-4;



%%% Environmental transmissibility params pulled from Breban
% (eta just inspired by Breban)
eta = [3 5 4 2]; % Decay rate spring summer fall winter
alpha = 1; % 'Environmental infectiousness' % Maybe just here for units?
rho = 10e-4; % Exposure rate
rho_w = rho; rho_d = rho; rho_h = rho;
omega = 10e5; % virus shedding rate
omega_w = omega; omega_d = omega; omega_t = omega;
eta_w = eta; eta_d = eta; eta_t = eta;
dose = 10e5; % number of virions for infectious dose

  % count when j==0;
% Number of sampe paths, time step, and end time
sim = 1000; outbreak=100;
sea=4;

% CTMC Model (Gillespie's algorithm)

for g=1:sea 
TI=90*(g - 1) + 60;
totext(g)=0;
clear x

for k=1:sim

% Initial Conditions


t=0; i_w =1; i_d =0; i_h =0; t_d =0; 
s_w =999; s_d =1000; s_h =1000;
r_w =0; r_d =0; r_h =0;
v_w = 0; v_d = 0; v_t = 0;
j=1;
N_w = 1000; N_d = 1000; N_h = 1000; 

% While loop
while  i_h +i_d +i_w +t_d  > 0 && i_h +i_d +i_w +t_d  < outbreak && N_w > 1 && N_d > 1 && N_h > 1 
%while j < 100000
  u2=rand; % Only need one random variable/ u2=rand; % two uniform random numbers

  beta_ww  = amp_bww*sin((t  + TI)*3.14/180 - 3.14) + beta_ww0;
  b_w = (0.5*tanh(57*sin(((t  + TI)*pi/180)-106*2*pi/360)-40)+0.5)/270;

  N_w=s_w +i_w +r_w ; % set population sizes
  N_d=s_d +i_d +r_d +t_d ;
  N_h=s_h +i_h +r_h ;

  ev1=b_w*dt; % birth of wild bird
  ev2=m_w*s_w *dt+ev1;% death of sw
  ev3=m_w*i_w *dt+ev2;% death of Iw
  ev4=dt*(beta_ww *s_w *i_w /N_w + rho_w*s_w *v_w /dose)+ev3; % infection of wild bird
  ev5=gamma_w*i_w *dt+ev4;% recovery of wild bird
  ev6=m_w*r_w *dt+ev5;% death of rw
  ev7=b_d*dt+ev6; %birth of chicken
  ev8=((beta_dd*s_d *i_d +beta_wd*s_d *i_w )/N_d + s_d *(rho_d*v_d +rho_d*v_w )/dose)*dt+ev7;% infection of sd by id or iw
  ev9=(beta_td*s_d *t_d /N_d + rho_d*s_d *v_t /dose)*dt+ev8;% infection of sd by td
  ev10=gamma_d*i_d *dt+ev9; % recovery of id
  ev11=gamma_d*t_d *dt+ev10;% recovery of td
  ev12=m_d*s_d *dt+ev11;% death of sd
  ev13=m_d*i_d *dt+ev12; % death of id
  ev14=m_d*t_d *dt+ev13;% death of td
  ev15=m_d*r_d *dt+ev14;% death of rd
  ev16=mu_dt*i_d *dt+ev15; % mutation of id
  ev17=b_h*dt+ev16;% birth of human
  ev18=((beta_th*s_h *t_d +beta_hh*s_h *i_h )/N_h + rho_h*s_h *v_t /dose)*dt+ev17;% infection of human
  ev19=m_h*s_h *dt+ev18; %death of sh
  ev20=m_h*i_h *dt+ev19;% death of ih
  ev21=m_h*r_h *dt+ev20;% death of rh
  ev22=gamma_h*i_h *dt+ev21; % recovery of human
  ev23=omega_w*i_w *dt/dose+ev22; % shed infectious dose v_w
  ev24=omega_d*i_d *dt/dose+ev23; % shed infectious dose v_d
  ev25=omega_t*t_d *dt/dose+ev24; % shed infectious dose v_t
  ev26=eta_w(g)*v_w *dt/dose+ev25; % decay infectious dose v_w
  ev27=eta_d(g)*v_d *dt/dose+ev26; % decay infectious dose v_d
  ev28=eta_t(g)*v_t *dt/dose+ev27; % decay infectious dose v_t
  ev29=1; % nothing happens
   
  t =t +dt; % interevent time

  if u2 <= ev1 % birth of wild bird
    %disp('birth of wild bird')
    s_w =s_w +1;

  elseif u2<=ev2  %death of sw
    %disp('death of sw')
    s_w =s_w -1;

  elseif u2<=ev3 %death of Iw
    %disp('death of Iw')
    i_w =i_w -1;

  elseif u2<=ev4  %infection of wild bird
    %disp('infection of wild bird')
    s_w =s_w -1;
    i_w =i_w +1;

  elseif u2<=ev5  %recovery of wild bird
    %disp('recovery of wild bird')
    i_w =i_w -1;
    r_w =r_w +1;

  elseif u2<=ev6  %death of rw
    %disp('death of rw')
    r_w =r_w -1;

  elseif u2<=ev7  %birth of chicken
    %disp('birth of chicken')
    s_d =s_d +1;

  elseif u2<=ev8  %infection of sd by id or iw
    %disp('infection of sd by id or iw')
    s_d =s_d -1;
    i_d =i_d +1;

  elseif u2<=ev9  %infection of sd by td
    %disp('infection of sd by td')
    s_d =s_d -1;
    t_d =t_d +1;

  elseif u2<=ev10  %recovery of id
    %disp('recovery of id')
    i_d =i_d -1;
    r_d =r_d +1;

  elseif u2<=ev11  %recovery of td
    %disp('recovery of td')
    t_d =t_d -1;
    r_d =r_d +1;

  elseif u2<=ev12 %death of sd
    %disp('death of sd')
    s_d =s_d -1;

  elseif u2<=ev13  %death of id
    %disp('death of id')
    i_d =i_d -1;

  elseif u2<=ev14  %death of td
    %disp('death of td')
    t_d =t_d -1;

  elseif u2<=ev15  %death of rd
    %disp('death of rd')
    r_d =r_d -1;

  elseif u2<=ev16  %mutation of id
    %disp('mutation of id')
    i_d =i_d -1;
    t_d =t_d +1;

  elseif u2<=ev17  %birth of human
    %disp('birth of human')
    s_h =s_h +1;

  elseif u2<=ev18  %infection of human
    %disp('infection of human')
    s_h =s_h -1;
    i_h =i_h +1;

  elseif u2<=ev19   %death of sh
    %disp('death of sh')
    s_h =s_h -1;

  elseif u2<=ev20 %death of ih
    %disp('death of ih')
    i_h =i_h -1;

  elseif u2<=ev21  %death of rh
    %disp('death of rh')
    r_h =r_h -1;

  elseif u2<=ev22  %Recovery of Ih
    %disp('recovery of ih')
    i_h =i_h -1;
    r_h =r_h +1;

  elseif u2<=ev23 % shed infectious dose v_w
    v_w =v_w +1;

  elseif u2<=ev24 % shed infectious dose v_d
    v_d =v_d +1;

  elseif u2<=ev25 % shed infectious dose v_t
    v_t =v_t +1;

  elseif u2<=ev26 % decay infectious dose v_w
    v_w =v_w -1;

  elseif u2<=ev27 % decay infectious dose v_d
    v_d =v_w -1;

  elseif u2<=ev28 % decay infectious dose v_t
    v_t =v_w -1;
    
 
  end
  j=j+1;
  
end  

end

if i_h +i_w +i_d +t_d ==0
    totext(g)=totext(g)+1;
end
end


sim_prob_ext_summer=totext(1)/sim
sim_prob_ext_fall=totext(2)/sim
sim_prob_ext_winter=totext(3)/sim
sim_prob_ext_spring=totext(4)/sim

R0h=beta_hh/(gamma_h+m_h)
R0w=beta_ww0/(gamma_w+m_w)
R0t=beta_td/(gamma_d+m_d)
R0d=beta_dd/(gamma_d+m_d+mu_dt)
toc
