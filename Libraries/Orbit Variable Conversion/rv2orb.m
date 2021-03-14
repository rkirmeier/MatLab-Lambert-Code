function [orb,Orbit] = rv2orb(R,V,u,Pram)
%% Position & Velocity Vector to Orbital Parameters
%   Detailed explanation goes here
... http://protean.io:3002/
%% Equations & Settings:
 ANG = @(A,B) acosd(dot(A,B)/(norm(A)*norm(B)));
             %orb = [a e i Omega w; vo]
%  Pram = {'a' 'e' 'i' 'O' 'w' 'vo'}; 
   orb = [0 0 0 0 0 0];
     Pram = {'a' 'e' 'i' 'O' 'w' 'vo'};
%% Standard:
 % Unit Vectot:
    I_hat = [1 0 0];
    J_hat = [0 1 0];
    K_hat = [0 0 1];
    
 % Magnitudes of R and V:
    r = norm(R);
    v = norm(V);
  
 % Specific Angular Momentum:
    H = cross(R,V);
      h = norm(H);     %h = absolute(H) = rp*vp = ra*va;

 % Node Vector:
    N = cross(K_hat,H);
      n = norm(N);    
  
%% Argument 1:
 % Ecentricity Vector:
    EV = (1/u) *((v^2 - u/r)*R - dot(R,V)*V);
      e = norm(EV);      %e = sqrt(1-(p/a));
      
%% Argument 2:
 % Semi-Latus Rectum:
    p = (h^2)/u;
 % Semi-major Axis:
    a = p/(1-e^2);   %u/((2*u/r)-v^2);
    b = sqrt((1-(e^2))*a^2);
     
%% Argument 3:
  % Inclination:                    (Angle between K and h)
     i = ANG(K_hat,H);  %acosd(dot(H,K_hat)/(h)); 
       if i>180
         i = i - 360;
       end     
     ...Always less than 180

%% Argument 4:
  % Longitude of the Ascending Node:  (Angle between I_hat and N)
     Omega = ANG(I_hat,N); %?  
       if dot(N,J_hat)<0
         Omega = 360 - Omega;
       end

%% Argument 5:
  % Argument of Periapsis:(Angle between the Acending Node & Periapsis)
     w = ANG(N,EV);    
       if dot(EV,K_hat)<0
         w = 360 - w;
       end
  % Longitude of Periapsis: (Angle from I_hat to Periapsis)
      II = Omega - w;
      
%% Argument 6:    
  % True Anonaly @ Epoch: (Angle Between Periapsis and Satellite)
      Vo = ANG(EV,R);%acos(r*e);   %ang(E,R,e,r);
        if dot(R,V)<0
          Vo = 360 - Vo;
        end
  % Argument of Latitude @ Epoch: (Angle between Acending Node and R)   
      uo = ANG(N,R);% w - Vo;
        if dot(R,K_hat)<0
          uo = 360 - uo;
        end
  % True Longitude @ Epoch: (Angle Between I_Hat and R)
      Lo = Omega + w + Vo;   %=II + Vo = Omega + uo
      
  % Period:
      TP = 2* pi/sqrt(u) * a^(3/2);
  % Time of Periapsis Passage: (Last time the satilite was at periapsis)[sec]
      T = -(TP/(2*pi)) * ( (Vo*pi/180) - atan( ((b-a)*sind(2*Vo)) / (b+a+(b-a)*cosd(2*Vo)) ) );
      
%% OTHER: 
  % Apo & Per:
      ra = p/(1-e); %a*(1+e)
      rp = p/(1+e); %a*(1-e)
    
  % Eccentric Anomaly:
      E = acos((e+cosd(Vo))/(1+e*cosd(Vo))); 

  % Mean Motion: [rad/s] 
      MM = sqrt(u/a^3);         
    
    
  % Mean Anomaly:
      M = E - e*sin(E);
  
%     M_deg = orb(6);
%       M_rad = deg2rad(M_deg);
%       E = M_rad; 
%       dE = 10^10;
%       eps = 1e-6; % [rad] control precision of Newton's method solution
%       while (abs(dE) > eps)
%         dE = (E - e * sin(E) - M_rad)/(1 - e * cos(E));
%         E = E -  dE;
%       end
%       
    
%% Output:
% switch count
%     case 0
        Orbit.R = R;    Orbit.r = r;
        Orbit.V = V;    Orbit.v = v;
        Orbit.H = H;    Orbit.h = h;
        Orbit.N = N;    Orbit.n = n;
        Orbit.EV = EV;  Orbit.e = e;
    
        Orbit.a = a;    Orbit.p = p;
        Orbit.i = i;
        Orbit.w = w;    Orbit.II = II;
        Orbit.Vo = Vo;  Orbit.uo = uo;
        Orbit.Lo = Lo;  Orbit.T = T;
        Orbit.ra = ra;  Orbit.rp = rp;
        Orbit.E = E;    Orbit.M = M;
        Orbit.TP = TP;  Orbit.MM = MM;
        Orbit.Omega = Omega;



  % Argument 1:
      switch Pram{(1)}
          case {'','a'}
             orb(1) = a;
          case 'p'
             orb(1) = p;
          otherwise
             warning('Wrong Value for Argument 1')
             orb(1) = 'err';
            return
      end
       
  % Argument 2:
      switch Pram{(2)}
          case {'','e'}
             orb(2) = e;
          otherwise
             warning('Wrong Value for Argument 2')
             orb(2) = 'err';
             return
      end
       
  % Argument 3:
      switch Pram{(3)}
          case {'','i'}
             orb(3) = i;
          otherwise
             warning('Wrong Value for Argument 3')
             orb(3) = 'err';
             return
      end
       
  % Argument 4:
      switch Pram{(4)}
          case {'','O','Omega','omega'}
             orb(4) = Omega;
          otherwise
             warning('Wrong Value for Argument 4')
             orb(4) = 'err';
             return
       end
           
  % Argument 5:
      switch Pram{(5)}
          case {'','o','w'}
             orb(5) = w;
          case {'II','PI','LOP'}
             orb(5) = II;
          otherwise
             warning('Wrong Value for Argument 5')
             orb(5) = 'err';
             return
       end
           
  % Argument 6:
      switch Pram{(6)}
          case {'','Vo','V0','vo','v0','TAE'}
             orb(6) = Vo;
          case {'uo','ALE'}
             orb(6) = uo;
          case {'lo','L','TLE','L'}
             orb(6) = Lo;
          otherwise
             warning('Wrong Value for Argument 6')
             orb(6) = 'err';
             return
      end
     
end

