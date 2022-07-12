%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Royce CTMC Model
% %%%%%%%%%%%%%%%%%%%%%%%%

function prob_ext_vec = OnePatchFast

N0=1000; %Initial population

% Parameters

amp_bww=-0.5; beta_ww0 = 0;
beta_wd=.51; beta_dd=.89; beta_td=.5; beta_th=.207; beta_hh=.078; %transmission rates
gamma_w=0.1444; gamma_d=0.1444; gamma_h=0.1444; %recovery rates
mu_dt=.499; %mutation to human transmissible 
m_w=.00082; m_d=.00082; m_h=.009; %death rates
%%% Environmental transmissibility params pulled from Breban
% (eta just inspired by Breban)
eta = [5 2 1 3]; % Decay rate summer, fall, winter, spring
alpha = 1; % "Environmental infectiousness" % Maybe just here for units?
rho = 0; % Exposure rate
omega = 10e5; % virus shedding rate

%%%
b_d=.0082*N0; b_h=N0*0.0118; %birth rates
dt=.00023844;

t = [0:dt:60];
totext_sum=0; totext_fall = 0; 
totext_spring = 0; totext_winter = 0;
% Number of sampe paths, time step, and end time
sim = 100; outbreak=200;
sea=4;

% CTMC Model (Gillespie's algorithm)
for g=1:sea 
    TI=91.25*g;

    for k=1:sim

        % Initialize vectors 
        
        j = 1;
        i_w=1;i_d=0; i_h=0; t_d=0; %initial infected
        r_w=0; r_d=0; r_h=0; %initial recovered
        N_w = N0; N_d = N0; N_h = N0; % because while loop conditions
        s_w=N_w-i_w; s_d=N_d-i_d-t_d; s_h=N_h-i_h; %initial susceptible
        x = 0;
        infected = 1; % for while loop condition
        

        % While loop
        %while  infected > 0 && infected < outbreak && N_w > 1 && N_d > 1 && N_h > 1 && j < 250000
        while j<100
            u1=rand; 
            func_val = sin((t(j)+TI)*(pi/180));
  
            beta_ww = amp_bww*func_val + beta_ww0;
            if func_val>0.5
                b_w=2.5;
            else
                b_w=0;
            end

             N_w=s_w+i_w+r_w; % set population sizes
             N_d=s_d+i_d+r_d+t_d;
             N_h=s_h+i_h+r_h;

             ev1=b_w*dt; % birth of wild bird
             ev2=m_w*s_w*dt+ev1;% death of sw
             ev3=m_w*i_w*dt+ev2;% death of Iw
             ev4=dt*s_w*(beta_ww*i_w/N_w + rho*(1 - exp(- alpha*omega*x)))+ev3; % infection of wild bird
             ev5=gamma_w*i_w*dt+ev4;% recovery of wild bird
             ev6=m_w*r_w*dt+ev5;% death of rw
             ev7=b_d*dt+ev6; %birth of chicken
             ev8=s_d*(beta_dd*i_d+beta_wd*i_w)/N_d*dt+ev7;% infection of sd by id or iw
             ev9=beta_td*s_d*t_d/N_d*dt+ev8;% infection of sd by td
             ev10=gamma_d*i_d*dt+ev9; % recovery of id
             ev11=gamma_d*t_d*dt+ev10;% recovery of td
             ev12=m_d*s_d*dt+ev11;% death of sd
             ev13=m_d*i_d*dt+ev12; % death of id
             ev14=m_d*t_d*dt+ev13;% death of td
             ev15=m_d*r_d*dt+ev14;% death of rd
             ev16=mu_dt*i_d*dt+ev15; % mutation of id
             ev17=b_h*dt+ev16;% birth of human
             ev18=s_h*(beta_th*t_d+beta_hh*i_h)/N_h*dt+ev17;% infection of human
             ev19=m_h*s_h*dt+ev18; %death of sh
             ev20=m_h*i_h*dt+ev19;% death of ih
             ev21=m_h*r_h*dt+ev20;% death of rh
             ev22=gamma_h*i_h*dt+ev21;

             if u1 <= ev1 % birth of wild bird
                 disp('birth s_w')
                 s_w=s_w+1;
    
             elseif u1<=ev2 %death of sw
                 disp('death s_w')
                 s_w=s_w-1;
    
             elseif u1<=ev3 %death of Iw
                 disp('death i_w')
                 i_w=i_w-1;
    
             elseif u1<=ev4 %infection of wild bird
                 disp('infection s_w')
                 s_w=s_w-1;
                 i_w=i_w+1;
    
             elseif u1<=ev5 %recovery of wild bird
                 disp('recovery i_w')
                 i_w=i_w-1;
                 r_w=r_w+1;
             elseif u1<=ev6 %death of rw
                 disp('death r_w')
                 r_w=r_w-1;
 
             elseif u1<=ev7 %birth of chicken
                 disp('birth s_d')
                 s_d=s_d+1;
    
             elseif u1<=ev8 %infection of sd by id or iw
                 disp('infection s_d')
                 i_d=i_d+1;
                 s_d=s_d-1;
    
             elseif u1<=ev9 %infection of sd by td
                 disp('infection s_d to t_d')
                 t_d=t_d+1;
                 s_d=s_d-1;
    
             elseif u1<=ev10 %recovery of id
                 disp('recovery i_d')
                 i_d=i_d-1;
                 r_d=r_d+1;
    
             elseif u1<=ev11 %recovery of td
                 disp('recovery t_d')
                 t_d=t_d-1;
                 r_d=r_d+1;

             elseif u1<=ev12 %death of sd
                 disp('death of s_d')
                 s_d=s_d-1;

             elseif u1<=ev13 %death of id
                 disp('death of i_d')
                 i_d=i_d-1;

             elseif u1<=ev14 %death of td
                 disp('death of t_d')
                 t_d=t_d-1;

             elseif u1<=ev15 %death of rd
                 disp('death of r_d')
                 r_d=r_d-1;

             elseif u1<=ev16 %mutation of id
                 disp('mutation of i_d')
                 i_d=i_d-1;
                 t_d=t_d+1;

             elseif u1<=ev17 %birth of human
                 disp('birth of human')
                 s_h=s_h+1;
 
             elseif u1<=ev18 %infection of human
                 disp('infection of s_h')
                 i_h=i_h+1;
                 s_h=s_h-1;
 
             elseif u1<=ev19 %death of sh
                 disp('death of s_h')
                 s_h=s_h-1;
 
             elseif u1<=ev20 %death of ih
                 disp('death of i_h')
                 i_h=i_h-1;
 
             elseif u1<=ev21 %death of rh
                 disp('death of r_h')
                 r_h=r_h-1;
 
             elseif u1<=ev22 %Recovery of Ih
                 disp('recovery of i_h')
                 i_h=i_h-1;
                 r_h=r_h+1;
 
             end
             val = exp(-eta(g)*dt);
             x = x*val + i_w*(1 - val)/eta(g);
             j=j+1;
            
            
            infected = i_h + i_w + i_d + t_d;
            if infected==0 && g==1
                %totext_sum=totext_sum+1;
            elseif infected==0 && g==2
                %totext_fall=totext_fall+1;
            elseif infected==0 && g==3
                %totext_winter=totext_winter+1;
            elseif infected==0 && g==4
                %totext_spring=totext_spring+1;
            end


        
        end  % end while loop

    end % end for loop for simulations
end % end for loop for seasons

Prob_ext_summer=totext_sum/sim;
Prob_ext_fall=totext_fall/sim;
Prob_ext_winter=totext_winter/sim;
Prob_ext_spring=totext_spring/sim;

prob_ext_vec = [Prob_ext_summer Prob_ext_fall Prob_ext_winter Prob_ext_spring];

end