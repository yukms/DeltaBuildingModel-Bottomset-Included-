%close all, clc
%% Getting input parameters first:
tic
%% Compute Mud Fall Velocity
% This subroutine implements the Dietrich relation for fall velocity.

Rep = ((Rm * g * Dm) ^ 0.5) * Dm / nu;
xx = log(Rep ^ 2) / log(10);
yy = -3.76715 + 1.92944 * xx - 0.09815 * (xx ^ 2) - 0.00575 * (xx ^ 3) + 0.00056 * (xx ^ 4);
Rf = ((10 ^ yy) / Rep) ^ (1 / 3);
vsm =  Rf * ((Rm * g * Dm) ^ 0.5);
        
% %% Compute Normal Flow
% % This subroutine computes the Shields stress, slope and depth of the normal flow
% % of the fluvial region associated with the volume sand feed rate per unit wifth qso,
% % the water discharge per unit width qw and the sand grain size Ds.
% 
% tsn = (qso / ((Rs * g * Ds) ^ 0.5 * Ds * al)) ^ (1 / nl) + tsc;
% Sn = Cza * (Rs * Ds * tsn) ^ 1.5 * g ^ 0.5 / qw;
% Hn = Rs * Ds * tsn / Sn;


%% Set Initial bed
%This subroutine sets up the initial bed so that the fluvial reach has 
% initial length rsi and slope Sfbi, and the subaqueous reach has initial 
% length rmax-rbi and slope Stbi.

for i = 1:1:N + 1
    sfbar(i) = dsbar * (i - 1);
    rfluv(i) = rsi * sfbar(i);
    etafluv(i) = etatI + rsi * Sfbi - rsi * sfbar(i) * Sfbi;
    r(i) = rfluv(i);
    eta(i) = etafluv(i);
end
rbi = rsi + (etatI - etabI) / Sa;
for i = 1:1:Ns + 1
    stbar(i) = dstbar * (i - 1);
    rturb(i) = rbi + (rmax - rbi) * stbar(i);
    etaturb(i) = etabI - Stbi * (rturb(i) - rbi);
    jj = N + 1 + i;
    r(jj) = rturb(i);
    eta(jj) = etaturb(i);
end
        etaturbi=etaturb;
rss = rsi;
rbb = rbi;
etatop = etatI;
etabot = etabI;
% etasdot = 0;
% etabdot = 0;

Hturb=Hturb-etaturb; %actually not turbidity zone, settling zone
Uturb = qw./Hturb;

time = 0;
ntrack = 0;
j = 0; 
% Finished = False;
% Bombed = False;
ystrata(:,j+1)=eta(:,1);
xstrata(:,j+1)=r(:,1);
shoreline(j+1,1)=rss;
etashoreline(j+1,1)=etatop;
fbb(j+1,1)=rbb;
etafbb(j+1,1)=etabot;
timeplot(j+1,1)=time;
%%

for j=1:1:Nprint
    %cibar=zeros(N+1,1);   %when to reset cibar to zeros ( j loop) - because this model runs for decades ( levee model ran for short period )
    
    %Inner loop
    for w = 1:1:Ntoprint
        %% Do Fluvial Backwater
        %This subroutine implements a numerical solution of the backwater 
        % equation on both the sand-bed and gravel-bed reaches

        Hfluv(N + 1) = xil - etafluv(N + 1);
        for i = 1:1:N
            fr2p = qw ^ 2 / g / Hfluv(N + 2 - i) ^ 3;  % Fr^2
            fnp = (etafluv(N + 1 - i) - etafluv(N + 2 - i) - Cfa * fr2p * rss * dsbar) / (1 - fr2p); %dH/dx^
            Hpred = Hfluv(N + 2 - i) - fnp; %  H predictor corrector
            fr2 = qw ^ 2 / g / Hpred ^ 3; %?
            fn = (etafluv(N + 1 - i) - etafluv(N + 2 - i) - Cfa * fr2 * rss * dsbar) / (1 - fr2);
            Hfluv(N + 1 - i) = Hfluv(N + 2 - i) - 0.5 * (fnp + fn); % Eueler step used 
            if Hfluv(N + 1 - i) < Hn
                Hfluv(N + 1 - i) = Hn; % Normal flow : Hn
            end
        end
        xil = xil + xiddot * dt;
        
        %% Find SandLoad and rsdot        %rsdot = S_s dot,shoreline moving rate (52)
        % This subroutine determines the sand transport rate qs(i) at every node in accordance with
        % (1), (5) and (12) of the notes, and also gets a first estimate of the migration speed
        % of the topset-foreset interface in accordance with (52).
        for i = 1:1:N+1
            ts = Cfa * (qw / Hfluv(i)) ^ 2 / (Rs * g * Ds);
            if ts <= tsc 
                qs(i) = 0;
            else
                qs(i) = (Rs * g * Ds) ^ 0.5 * Ds * al * (ts - tsc) ^ nl;
            end
            
            
%             if  i==1
%                  qms(i)=qmo;
%              else
%                  if i<N+1
%                      qms(i)=qms(i-1) - Lamms*(qs(i-1)-qs(i));
%                      if qms(i)<0
%                          qms(i)=0;
%                      end
%                  end
%              end  %% mud deposited per sand  in topset 
             
             
        end
        rsdot(j,1) = 1 / Sa *intM* (qs(N + 1)*(1+Lamms) / ((1 - lps) * (rbb - rss))); %S_s dot (52) :shock condition
%         qms(N+1) =qms(N) - qs(N+1) * Lamms;  %% mud deposited per sand in foreset
%         if qms(N)<0
%             qms(N)=0;
%         end
        

        
        if (j == 1) && (w == 1)
            rbdot(j,1) = rsdot(j,1) * Sa / (Sa - Sturb(1));
        else
            rbdot(j,1) = rsdot(j,1) + 1 / Sa * (etasdot(j,1) - etabdot(j,1));
        end 
        


   

        %% Find New etafluv and rss
        % This subroutine implements the Exner equation (47) of the notes over the fluvial
        % region to find the new bed elevation.  In addition, it updates the migration speed
        % and position of the topset-foreset interface according to (52) by adding a term.
        etatopold = etatop;
        for i = 1:1: N + 1
            if i == 1 
                qsdev(i,1) = (qs(i) - qso) / dsbar;
                etafdev(i,1) = (etafluv(i + 1) - etafluv(i)) / dsbar;
            else
                if i == N + 1
                    qsdev(i,1) = (qs(i) - qs(i - 1)) / dsbar;
                    etafdev(i,1) = (etafluv(i) - etafluv(i - 1)) / dsbar;
                else
                    qsdev(i,1) = (qs(i) - qs(i - 1)) / dsbar;
                    etafdev(i,1) = (etafluv(i + 1) - etafluv(i)) / dsbar;
                end
            end
            etafluv(i) = etafluv(i) + intM*(dt * (-qsdev(i,1) / rss / (1 - lps)*(1+Lamms))) + (sfbar(i) * rsdot(j,1) / rss * etafdev(i,1));
            eta(i) = etafluv(i);  %Exner equation
        end
        etatop = etafluv(N + 1);
        etasdot(j,1) = (etatop - etatopold) / dt;
        rsdot(j,1) = rsdot(j,1) - etasdot(j,1) / Sa;
        rss = rss + rsdot(j,1) * dt;
        for i = 1:1: N + 1
            rfluv(i) = sfbar(i) * rss;
            r(i) = rfluv(i);
        end

             
        
        
        %% Find New etaturb and rbb
        % This subroutine implements the Exner equation (51) of the notes over the turbidity
        % current region to find the new bed elevation.  It also updates the migration speed
        % and position of the foreset-bottomset interface with the aid of (53) of the notes.
        etabotold = etabot;
%         cibar=zeros(N+1,1);   %when to reset cibar to zeros ( j loop) - because this model runs for decades ( levee model ran for short period )
     complex2=1/(rmax-rbb);

        for i = 1:1: Ns + 1
            jj = N + 1 + i;
            
            if i == Ns + 1
                etatdev = (etaturb(i) - etaturb(i - 1)) / dstbar;
            else
                if i == 1
                    etatdev = (etaturb(i + 1) - etaturb(i)) / dstbar;
                else
                    etatdev = (etaturb(i + 1) - etaturb(i)) / (dstbar);
                end
            end
%                             qm(1)=qms(N+1);
                              qm_ghost=qmo;
                            complex(i)=rbdot(j,1)*(1-stbar(i))/(rmax-rbb);
        end
        
        for i = 1:1:Ns+1
            jj=N+1+i;
            
            if i == 1 % dt :  0.00025 * timeyr  ,  too large ?  when dt : 0.00025 is better
                            dqm(i)= cibar(i)*Uturb(i)*Hturb(i)-qm_ghost;
%                             dqm(i)= (cibar(i)*Uturb(i)*Hturb(i)-qmo)*dsbar;
%                             cibar(i)=cibar(i) +( complex(i) *(dqm(i)/qw)/dstbar + ( -vsm * ro * cibar(i) /(1-lpm) - complex2 * dqm(i)/dstbar ) / Hturb(i) )*dt; % without H change _ dqm used instead of cibar dev
                            cibar(i) = cibar(i) + ( complex(i)*Hturb(i)*(dqm(i)/qw)/dstbar - complex2*dqm(i)/dstbar - vsm*ro*cibar(i) -complex(i)*cibar(i)*etatdev +cibar(i)*deta(i)/dt )*dt / Hturb(i); %changing dH
            else
%                             dqm(i)=((cibar(i-1)-cibar(i))*Uturb(i)*Hturb(i))*dsbar;
                            dqm(i)=(cibar(i)-cibar(i-1))*Uturb(i)*Hturb(i);
                            
%                             cibar(i)=cibar(i) +( complex(i) *(cibar(i)-cibar(i-1))/dstbar + ( -vsm * ro * cibar(i)/(1-lpm) - complex2 * dqm(i)/dstbar ) / Hturb(i) )*dt; % without H change _ dqm used instead of cibar dev
                            cibar(i) = cibar(i) + ( complex(i)*Hturb(i)*(cibar(i)-cibar(i-1))/dstbar - complex2*dqm(i)/dstbar - vsm*ro*cibar(i) -complex(i)*cibar(i)*etatdev + cibar(i)*deta(i)/dt )*dt / Hturb(i); % changing dH

            end
           
%             qm=Uturb.*Hturb.*cibar;
        end
        
        
        for i=1:1:Ns+1
            jj=N+1+i;
%                  if i==1
%                 qm(i)=qms(N+1)+dqm(i);
%               else
%                 qm(i) = qm(i-1)+dqm(i);
%               end  % qm is higher with this code
%             
            qm(i)=cibar(i)*Uturb(i)*Hturb(i);
              
%                   if ro==0
%               etadeb=0;
%             end
%               
              
            etaturb(i) = etaturb(i) + rbdot(j,1) * (1 - stbar(i)) / (rmax - rbb) * etatdev * dt; %term 2
            deta(i)=(rbdot(j,1) * (1 - stbar(i)) / (rmax - rbb) * etatdev * dt);
            etaturb(i) = etaturb(i) + 1 / (1 - lpm) * vsm * ro * qm(i) / Uturb(i) / Hturb(i) * dt; %term 1 
            deta(i) = deta(i)  + (1 / (1 - lpm) * vsm * ro * qm(i) / Uturb(i) / Hturb(i) * dt); 
      
        
            
            eta(jj) = etaturb(i);
            Hturb(i) = Hturb(i)-deta(i);
            Uturb(i)=qw/Hturb(i);
        end
%             qm(i+1) = qm(i) -(rmax-rbb)*vsm*ro*qm(i)/(Uturb(i)*Hturb(i));
%             qm(i+1)=(qm(i)-qw*(deta(i)-rbdot*(1-stbar(i))/(rmax-rbb)*etatdev)
%             qm(i+1)=qm(i)-ro*vsm*qm(i)/(Hturb(i)*Uturb(i))*dsbar;  %is this right ?  &dsbar?
%             qm(i+1)=qm(i)*0.9;
%             qm(i+1)=qm(i) - (deta(i)*(1-lpm)/vsm/ro)*qw;
%         end
        
        
        
        etabot = etaturb(1);
        etabdot(j,1) = (etabot - etabotold) / dt;   % always be positive in Purely depositional Sys
        rbdot(j,1) = rsdot(j,1) + 1 / Sa * (etasdot(j,1) - etabdot(j,1));
        rbb = rbb + rbdot(j,1) * dt;
        if rbb > rmax
            disp('Toe prograded past rmax');
            break;
            
        else
            for i = 1:1: Ns + 1
                jj = N + 1 + i;
                rturb(i) = rbb + (rmax - rbb) * stbar(i);
                r(jj) = rturb(i);
            end
            time = time + dt;
            ntrack = ntrack + 1;
        end
        
%                      washload=washload + qm(Ns+1)*dt; 
%                      capturemud=capturemud + (qmo-qms(N+1))*dt;

                     
                     
    end 
    ystrata(:,j+1)=eta(:,1);
    xstrata(:,j+1)=r(:,1);
    shoreline(j+1,1) = rss;
    etashoreline(j+1,1) = etatop;
    fbb(j+1,1) = rbb;
    etafbb(j+1,1) = etabot;
    timeplot(j+1,1) = time/timeyr;
%     timeplot(j+1,1)= time/timeyear;
end

%% Plot output
figure, plot(xstrata, ystrata, 'k');
hold; plot(rfluv, Hfluv+etafluv, 'b');
titletext1 = ['Strata ' num2str(time/timeyr) ' yr-run ' num2str(xiddot*1000*timeyr) ' mm/yr-SLR'];
title(titletext1); xlabel('distance [m]'); ylabel('elevation [m]');


plot(shoreline, etashoreline, '-o');
plot(fbb, etafbb, '-v');
ylim([0 21])
saveas(gcf,figurename2,'jpg')

% % time versus shoreline location
% 
% titletext2 = ['Shoreline and Foreset-Bottomset Trajectories'];
% title(titletext2); xlabel('distance [m]'); ylabel('elevation [m]');
% 
% figure, plot(timeplot, shoreline,  '-o');
% hold; plot(timeplot, fbb, '-v');
% titletext2 = ['Shoreline and Foreset-Bottomset Trajectories'];
% title(titletext2); xlabel('time [yr]'); ylabel('distance [m]');
% toc
% 
% % 
% % 
% 
% % Dimensionless shoreline fbb figure
% 
% Lscale = ((qso+qmo)) / (xiddot) ; % m 
%  Tscale = Lscale / (xiddot) /  (1/Sfbi - 1/Sa ) /timeyr ; %  yr
%  figure(2), plot(shoreline./Lscale, timeplot./Tscale,'-k');
%  hold; plot(fbb./Lscale, timeplot./Tscale,'--k');
%  titletext2 = ['Shoreline and Foreset-Bottomset Transition Trajectories'];
%  title(titletext2); ylabel('time [t*]'); xlabel('distance [L*]');
% saveas(gcf,figurename1,'jpg')

%% wash load cal
% washload=washload/((qmo*dt*Ntoprint*Nprint)-capturemud) ;
% ..   'qm' of last node/input-captured ...

%% topset foreset length cal
for l = 1 : Nprint
    topsetlength(l)=sqrt((xstrata(N+1,l)-xstrata(1,l))^2 + (ystrata(N+1,l)-ystrata(1,l))^2);
    foresetlength(l)=sqrt((xstrata(N+2,l)-xstrata(N+1,l))^2 + (ystrata(N+2,l)-ystrata(N+1,l))^2);
end
    
