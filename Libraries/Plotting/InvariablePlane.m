clear,clc
digits (16);
format long;
% AerE 551, fall 2020
% Ryan Kirmeier
% The solar system’s invariable plane
% aanda.org/articles/aa/pdf/2012/07/aa19011-12.pdf#page=5&zoom=100,0,-84
%% Input:
  T1 = datetime(2020,12,01);
  EPH = '405';
  Center = 'SolarSystem';
  N = 11;       %Number of Bodies
  
%% Constants:
  G = 6.674083131313131313e-11;     %Gravitational Constent [m^3/(kg?s^2)]
 %Unit Vector:
  ih = [1 0 0]; 
  jh = [0 1 0]; 
  kh = [0 0 1];
  
%Bodies:
% SUN:
    m(1) = 19890e30;     % [kg]
      rSun = 6.957e8;      % [m]  Sun Radius @ Equator
    [R(1,1:3),V(1,1:3)] = ...
      planetEphemeris(juliandate(T1),'Sun',Center,EPH,'km'); 
% MERCURY:
    m(2) = 0.33e30;      % [kg]
    [R(2,1:3),V(2,1:3)] = ...
      planetEphemeris(juliandate(T1),'MERCURY',Center,EPH,'km'); 
% VENUS:
    m(3) = 4.87e30;      % [kg]
    [R(3,1:3),V(3,1:3)] = ...
      planetEphemeris(juliandate(T1),'VENUS',Center,EPH,'km');
% EARTH:
    m(4) = 5.97e30;      % [kg]
     rEarth = 6378e3;      % [m]  Earth Radius @ Equator
    [R(4,1:3),V(4,1:3)] = ...
      planetEphemeris(juliandate(T1),'EARTH',Center,EPH,'km');
% MOON:
    m(5) = 0.073e30;     % [kg]
    [R(5,1:3),V(5,1:3)] = ...
      planetEphemeris(juliandate(T1),'Moon',Center,EPH,'km');
% MARS:
    m(6) = 0.642e30;     % [kg]
     rMars = 3396e3;       % [m]  Mars Radius @ Equator
    [R(6,1:3),V(6,1:3)] = ...
      planetEphemeris(juliandate(T1),'MARS',Center,EPH,'km');
% JUPITER:
    m(7) = 1898e30;      % [kg]
    [R(7,1:3),V(7,1:3)] = ...
      planetEphemeris(juliandate(T1),'JUPITER',Center,EPH,'km');
% SATURN:
    m(8) = 568e30;       % [kg]
    [R(8,1:3),V(8,1:3)] = ...
      planetEphemeris(juliandate(T1),'SATURN',Center,EPH,'km');
% URANUS:
    m(9) = 86.8e30;      % [kg]
    [R(9,1:3),V(9,1:3)] = ...
      planetEphemeris(juliandate(T1),'URANUS',Center,EPH,'km');
% NEPTUNE:
    m(10) = 102e30;      % [kg]
    [R(10,1:3),V(10,1:3)] = ...
      planetEphemeris(juliandate(T1),'NEPTUNE',Center,EPH,'km');
% PLUTO:
    m(11) = 0.0146e30;   % [kg]
    [R(11,1:3),V(11,1:3)] = ...
      planetEphemeris(juliandate(T1),'PLUTO',Center,EPH,'km');
										
R = R*10^3;
     
%% Equations:
    u = @(M) G*M;                        %µ
    fL = @(m,r,v) m * cross(r,v);
      fL1 = @(m,r,v) m * (r(2)*v(3)-r(3)*v(2));
      fL2 = @(m,r,v) m * (r(3)*v(1)-r(1)*v(3));
      fL3 = @(m,r,v) m * (r(1)*v(2)-r(2)*v(1));
%     Ltot = sqrt(L1^2 + L2^2 + L3^2);
%% INITIAL Values:
    mu = u(m(1));                %(m^3/s^2)
    L1 = 0; L2 = 0; L3 = 0;
    for j = 1:N
        L1 = L1 + fL1(m(j),R(j,1:3),V(j,1:3));
        L2 = L2 + fL2(m(j),R(j,1:3),V(j,1:3));
        L3 = L3 + fL3(m(j),R(j,1:3),V(j,1:3));
    end
        
    Ltot = sqrt(L1^2 + L2^2 + L3^2);
    
    i = acosd(L3/Ltot);
    O = atan(-L1/L2);
%     L1 = Ltot*sin(?)*sin(i);
%     L2 = -Ltot*cos(?)*sin(i);
%     L3 = Ltot*cos(i);
     
  ICRF_TRANS = TransMat(-i,O,0);
    
    
    
    
    
% % SUN:
%         mSUN = 19890e30;     % [kg]
% % MERCURY:
%     mMERCURY = 0.33e30;      % [kg]
% % VENUS:
%       mVENUS = 4.87e30;      % [kg]
% % EARTH:
%       mEARTH = 5.97e30;      % [kg]
% % MOON:
%        mMOON = 0.073e30;     % [kg]
% % MARS:
%        mMARS = 0.642e30;     % [kg]
% % JUPITER:
%     mJUPITER = 1898e30;      % [kg]
% % SATURN:
%      mSATURN = 568e30;       % [kg]
% % URANUS:
%      mURANUS = 86.8e30;      % [kg]
% % NEPTUNE:
%     mNEPTUNE = 102e30;       % [kg]
% % PLUTO:
%       mPLUTO = 0.0146e30;    % [kg]