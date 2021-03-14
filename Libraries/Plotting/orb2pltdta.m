function [Orbit,Position,Velocity] = orb2pltdta(DIM,SAT,rPlanet,mu,toll)
%   Orbit Plot Values
%% *History*
%  When       Who       What
%  ---------- -------- --------------------------------------------------
%  2019/07/11 mnoah     original code
%  2020/21/10 kirmeier  the butcherer of original code

%   Detailed explanation goes here
  
%% Input: orb = [a e i Omega w; vo]
  toll = 100;
    a = SAT.a;
    p = SAT.p;
    e = SAT.e;
    i = SAT.i;
    O = SAT.Omega;
    w = SAT.w;
    Vo = SAT.Vo;
    
%     E = SAT.E; %acos((e+cosd(orb(6)))/(1+e*cosd(orb(6))));%2 * atan(  tand(orb(6)/2)  *  sqrt((1+e)/(1-e))  );       %5.672438320452783;   %
       
%     n = SAT.MM; %sqrt(mu/a^3);         % [rad/s] mean motion
%%    
% %For Mean Anomaly:
%     M_deg = orb(6);
%       M_rad = deg2rad(M_deg);
%       E = M_rad; 
%       dE = 10^10;
%       eps = 1e-6; % [rad] control precision of Newton's method solution
%       while (abs(dE) > eps)
%         dE = (E - e * sin(E) - M_rad)/(1 - e * cos(E));
%         E = E -  dE;
%       end
%% *Plot the Plane of the Orbit*   
%     dMdt = n;                               %dM/dt rad/s
%     dEdt = dMdt/(1 - e*cos(E));             %dE/dt rad/s
%     dpdt = -a*sin(E)*dEdt;                  %dp/dt m/s [rad/s] p component velocity
%     dqdt = a*cos(E)*dEdt*sqrt(1 - e^2);     %dq/dt m/s [rad/s] q component velocity
%   %Plot Body: 
%     %plotbody(3,[0 0 0],rPlanet,0.01,[0.75 1.00 0.75]);
%     
%   %Plot Orbit:    
%    %plot(pvals,qvals,'.b');
%      Evals = 0:(1/toll):360.0; % [deg] values of the eccentric anomaly around orbit 
%        pval = a*(cosd(Evals)-e); % [m] orbit positions
%        qval = a*sqrt(1 - e^2)*sind(Evals); % [m] orbit positions
%          R_PQWa(:,1) = pvals.';
%          R_PQWa(:,2) = qvals.';
%          R_PQWa(:,3) = 0;

%   %Plot Position:    
%     %plot(p,q,'.r','MarkerSize',25);  
%       p = a*(cos(E) - e);          %[m] coordinate along axis through center and perigee
%       q = a*sqrt(1 - e^2)*sin(E);  %[m] coordinate passing through focus and perpendicular to p-axis
%         Position = [p,q];                    %[m] coordinate
%% 3D

    R_Trans = TransMat(i,O,w);
    
%   V_PQW = zeros(360*toll,3);
%   V_IJK = zeros(360*toll,3);
    R_PQWa = zeros(360*toll,3);
    R_PQWb = zeros(360*toll,3);
    R_IJK = zeros(360*toll,3);
  
    for anom = 1:1:(361*toll) % [deg] values of the True anomaly around orbit
        
      %R_PQWa:
       Rvalpa = a*(cosd(anom)-e); % [m] orbit positions
       Rvalqa = a*sqrt(1 - e^2)*sind(anom); % [m] orbit positions
         R_PQWa(anom,1) = Rvalpa.';
         R_PQWa(anom,2) = Rvalqa.';
         R_PQWa(anom,3) = 0;

      %R_PQWb:
        rval = p/(1+e*cosd(anom/toll));
          Rvalp = rval*cosd(anom/toll);
          Rvalq = rval*sind(anom/toll);
            R_PQWb(anom,1) = Rvalp;
            R_PQWb(anom,2) = Rvalq;
            R_PQWb(anom,3) = 0;
            
      %R_IJK = R_Trans*R_PQW;
          R_ijk = R_Trans.*(R_PQWa(anom,1:3));
            R_IJK(anom,1) = R_ijk(1,1)+R_ijk(1,2)+R_ijk(1,3);
            R_IJK(anom,2) = R_ijk(2,1)+R_ijk(2,2)+R_ijk(2,3);
            R_IJK(anom,3) = R_ijk(3,1)+R_ijk(3,2)+R_ijk(3,3);      
    end
    

%             
%       %V_IJK = R_Trans*V_PQW;
%           V_ijk = R_Trans*R_PQW(T_anom);
%             V_IJK(T_anom,1) = V_ijk(1);
%             V_IJK(T_anom,2) = V_ijk(2);
%             V_IJK(T_anom,3) = V_ijk(3);
    
    

        
%    %Plot Arrow:
%      %Plot Arrow(Line): 
%        %plot([p p1],[q q1],'r');
%           p1 = p+dpdt*300;
%           q1 = q+dqdt*300;
%             a1 = [p1,q1];
%         
%       %Plot Arrow(Head):
%          theta = -atan2d(dqdt,dpdt);   
%         %Arrow Head (Inside)
%            %plot(xl,yl,'r');
%              xi = p1 + [0 -0.1*rPlanet*cosd(theta+30)];
%              yi = q1 + [0  0.1*rPlanet*sind(theta+30)];
%         %Arrow Head (Outside)
%            %plot(xr,yr,'r');
%              xr = p1 + [0 -0.1*rPlanet*cosd(theta-30)];
%              yr = q1 + [0  0.1*rPlanet*sind(theta-30)];
%              
%                   ah = {xi,yi,xr,yr};
              
%% Output:
    switch DIM
        case 2
            
            Orbit = R_PQWb; 
            
            %V_PQW:
              vval = sqrt(mu/p);
                Vvalp = vval*-sin(Vo);
                Vvalq = vval*(e + cos(Vo));
                  V_PQW(1) = Vvalp;
                  V_PQW(1) = Vvalq;
                  V_PQW(1) = 0;
            Velocity = V_PQW; 
            
           %Position:              
                Rvalp = SAT.r*cos(Vo);
                Rvalq = SAT.r*sin(Vo);
                  R(1) = Rvalp;
                  R(2) = Rvalq;
                  R(3) = 0;
            Position = R;
            
        case 3
            Orbit = R_IJK;
            Velocity = SAT.V;  
            Position = SAT.R;  
    end
            


end


%     E_deg_epoch = rad2deg(E); 
%     n_deg_per_s = rad2deg(n); % [deg/s] mean motion

% *Rotate To ECI*
%    Rz_Omega = [ ...
%     [cosd(orb(4)) sind(orb(4)) 0]; ...
%     [-sind(orb(4)) cosd(orb(4)) 0]; ...
%     KH];
%    Rx_i = [ ...
%     IH; ...
%     [0 cosd(orb(3)) sind(orb(3))]; ...
%     [0 -sind(orb(3)) cosd(orb(3))]];
% 
%    Rz_omega = [ ...
%     [cosd(orb(5)) sind(orb(5)) 0]; ...
%     [-sind(orb(5)) cosd(orb(5)) 0]; ...
%     KH]; 
% 
%   [Year,Month,Day,H,M,S] = datevec(OE.epoch);
% 
%   [Year,Month,Day,H,M,S] = datevec(OE.epoch);
%   HourUTC = H + M/60.0 + S/3600.0;
%   jd = juliandate(Year,Month,Day,HourUTC,0,0);
%   jd0 = juliandate(Year,Month,Day,0,0,0); %Form time in Julian centuries from J2000
%   T = (jd - 2451545.0d0)./36525.0d0;