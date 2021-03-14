clear,clc,close all
digits (16);
format long;
% AerE 551, fall 2020
% Ryan Kirmeier
% Generation of The Pork-Chop Plot:
%% Input:
  Pram = {'a' 'e' 'i' 'Omega' 'w' 'Vo'};   %Orbital Element Choice for orb Function
  Method = 5;...'Universal Variable';        %1)Universal Variable Method 2)Lagrange Method 
  IttMethod = 'Bisection';
  tolerance = 1e-5;              %Tollerance
  Dim = 3;
  EPH = '432t';
  Center = 'SolarSystem';     
  radi = 410e3;
  %Departure Date
    t11 = datetime(2011,09,10,0,0,0);...2011,09,10,0,0,0...2007,07,27,0,0,0
    t12 = datetime(2012,01,08,0,0,0);...2012,01,08,0,0,0...2007,11,24,0,0,0
  %Arrival Date
    t21 = datetime(2012,04,20,0,0,0);...2012,04,20,0,0,0...2008,02,17,0,0,0
    t22 = datetime(2013,02,14,0,0,0);...2013,02,14,0,0,0...2009,02,11,0,0,0
  %Tollerance:
    step = 1;
%% Other Data:       
 %Bodies:
    dt = t11 +days(1);
   %Sun:
     Body_Center = bodysort();
       Body_Center.name = 'Sun';
       Body_Center.mass = (1.989e30);   % [kg] Sun Mass
       Body_Center.radius=6.957e8;      % [m]  Sun Radius @ Equator
   %Earth:
     Body_Start = bodysort();
       Body_Start.name = 'Earth';
       Body_Start.mass = 5.9722e24;     % [kg] Earth Mass
       Body_Start.radius=6378e3;        % [m]  Earth Radius @ Equator
   %Mars:
     Body_Target = bodysort();
       Body_Target.name = 'Mars';
       Body_Target.mass=0.64171e24;...   % [kg] Mars Mass
       Body_Target.radius=3396e3;...     % [m]  Mars Radius @ Equator
       
% %      mu = (6.674083131313131313e-11)*(Body_Center.mass); 
         mu = 132712441933.0;% [m^3/s^2]
  %Invariable Plane Transfer Matrix:
     ICR2_TRANS1 = ICRF2IVP(t11,11,EPH,Center);
     ICR2_TRANS2 = ICRF2IVP(dt,11,EPH,Center);
     
  %Body Info:
   fprintf('Kepler Elements for Planetary Bodies:[%s %s %s %s %s; %s]\n',...
       Pram{(1)},Pram{(2)},Pram{(3)},Pram{(4)},Pram{(5)},Pram{(6)});  
     for body = ['C' 'S' 'T']
      %Grab Data:
         switch body
           case 'C'
              Body = Body_Center;
           case 'S'
              Body = Body_Start;
              T = t11;
           case 'T'
              Body = Body_Target;
              T = t21;
         end
            
      %Get Body Info:
         [Body.P1,Body.V1] =...
            planetEphemeris(juliandate(t11),Center,Body.name,EPH,'km'); 
         [Body.P2,Body.V2] =...
            planetEphemeris(juliandate(dt),Center,Body.name,EPH,'km'); 
         
         Body.P1= matsolv(Body.P1,ICR2_TRANS1);...*10^3; 
         Body.V1= matsolv(Body.V1,ICR2_TRANS1);...*10^3;
         Body.P2= matsolv(Body.P2,ICR2_TRANS2);...*10^3; 
         Body.V2= matsolv(Body.V2,ICR2_TRANS2);...*10^3;
            
         [Body.E6, Body.Initial] = rv2orb(Body.P1, Body.V1, mu,Pram); 
         [Body.E6, Body.Final  ] = rv2orb(Body.P2, Body.V2, mu,Pram);
                  Body.Avg = AvgKepElm(Body.Initial, Body.Final);
                  Body.MM = Body.Initial.MM;
                  Body.M = Body.Initial.M;
         %Time:
                  Body.T = t11 + seconds(Body.Initial.T);
         
      %Plot Elements:
         [Body.Path1, Body.P11, Body.V11] =...
            orb2pltdta(Dim,Body.Initial,Body.radius,mu);
         [Body.Path2,Body.P21,Body.V21] =...
            orb2pltdta(Dim,Body.Final,  Body.radius,mu);
         [Body.Path3,~,~] =...
            orb2pltdta(Dim,Body.Avg,    Body.radius,mu);
         

       fprintf('  %s:\n',Body.name);  
         fprintf("    Initial Orbit = [%.5g,%.5f,%.3f\260,%.3f\260,%.3f\260;%.3f\260]\n",...
             Body.Initial.a,    Body.Initial.e,Body.Initial.i,...
             Body.Initial.Omega,Body.Initial.w,Body.Initial.Vo);
         
         fprintf("    Final Orbit   = [%.5g,%.5f,%.3f\260,%.3f\260,%.3f\260;%.3f\260]\n",...
             Body.Final.a,    Body.Final.e,Body.Final.i,...
             Body.Final.Omega,Body.Final.w,Body.Final.Vo);
                             
%          Body.delta_Vo = Body.Avg.Vo_2 - Body.Avg.Vo_1;
      %Save Data:
         switch body
           case 'C'
              Body_Center = Body;
           case 'S'
              Body_Start = Body;
           case 'T'
              Body_Target= Body;
         end
     end
        fprintf('\n\n');
clear ICR2_TRANS2 u G body Body
%% Iterations:
     radi = radi+Body_Start.radius;
     MMS= Body_Start.MM;
     MS = Body_Start.M;
     TS = Body_Start.T;
     e1 = Body_Start.Avg.e;
% %      p1 = Body_Start.Avg.p;
% %      vval1 = sqrt(mu/p1);
     
     MMT= Body_Target.MM;
     MT = Body_Target.M;
     TT = Body_Target.T;
     e2 = Body_Target.Avg.e;
% %      p2 = Body_Target.Avg.p;
% %      vval2 = sqrt(mu/p2);
     
% %     TransS = TransMat(Body_Start.Avg.i,...
% %         Body_Start.Avg.Omega,Body_Start.Avg.w);
% %     TransT = TransMat(Body_Target.Avg.i,...
% %         Body_Target.Avg.Omega,Body_Target.Avg.w);
% %     Vof = @(M,e) M...
% %                  +( (-1/4)*(e^3)+2*e *sin(1*M))...
% %                  +(  (5/4)*(e^2)     *sin(2*M))...
% %                  +((13/12)*(e^3)     *sin(3*M));
  %Matrix Size:
     [~,~,daysx,~,~,~]...
         = datevec(between(t11,t12,'Days'));  
     [~,~,daysy,~,~,~]...
         = datevec(between(t21,t22,'Days'));
     C3 = zeros(daysx/step,daysy/step);
     DECLd = C3;
     RASCd = C3;
     TOF = C3;
     DVT = C3;
     Vinf = C3;
  
  row = 0;
  runtime = C3;
  a=tic;
  
for T1 = t11:days(step):t12
   %Time:
      row = row + 1;
      column = 0;
      disp(row);
   %Position:              
      [RS,VS] =...
         planetEphemeris(juliandate(T1),Center,Body_Start.name,EPH,'km');
% %            RS=matsolv(RS,ICR2_TRANS1);... (10^3)*RS;... 
% %            VS=matsolv(VS,ICR2_TRANS1);... (10^3)*VS;...
           
    for T2 = t21:days(step):t22
       column = column + 1;
       b=tic;          
       %Position:              
          [RT,VT] =...
             planetEphemeris(juliandate(T2),Center,Body_Target.name,EPH,'km');
% %                RT = matsolv(RT,ICR2_TRANS1);... (10^3)*RS;... 
% %                VS = matsolv(VT,ICR2_TRANS1);... (10^3)*VS;...
       %Time:
          [~,~,tod,~,~,~]= datevec(between(t11,T1,'Days'));
                Tx(row,column) = tod;
          [~,~,tod,~,~,~]= datevec(between(t21,T2,'Days'));
                Ty(row,column) = tod;
          [~,~,~,toh,tom,tos] = datevec(between(T1,T2,'time'));        %Time of Flight in hours
          TOFs = 60*(60*toh + tom) + tos;   %Time of Flight in Seconds
          N = round((TOFs/Body_Target.Initial.TP)-.5);
          
       %Velocity:   
          [Vd,Va,~,~,~,~]= lambert(mu,RS,RT,TOFs,0,Method,IttMethod,tolerance,Body_Center);
        DVd(1) = (Vd(1) - VS(1));.../((10e3));
        DVd(2) = (Vd(2) - VS(2));.../((10e3));
        DVd(3) = (Vd(3) - VS(3));.../((10e3));
        
        DVa(1) = (VT(1) - Va(1));.../((10e3));
        DVa(2) = (VT(2) - Va(2));.../((10e3));
        DVa(3) = (VT(3) - Va(3));.../((10e3));
        
        C3(row,column) = (norm(DVd))*(norm(DVd));...;
        Vinf(row,column) = norm(DVa);
        TOF(row,column) = TOFs/((60^2)*24);
        DVT(row,column) = norm(DVd) + norm(DVa);
        
%         if (C3(row,column) <= 50)
%               DECLd(row,column) = 90 - acosd(DVd(3) / norm(DVd));
%               RASCd(row,column) = atan2d(DVd(2), DVd(1));
%         else
%               DECLd(row,column) = NaN;
%               RASCd(row,column) = NaN;
%         end
% %         runtime(row*(300)+column-1) = toc(b);        % compute arrival dla and rla
% %         toc(b);
    end
end
runtimeavg = toc(a);

  C3_lvl = [01,02,03,04,05,06,07,08,09,10 ...
           ,11,12,13,14,15,16,17,18,19,20 ...
           ,21,22,23,24,25,26,27,28,29,30 ...
            ];
figure(1)
 [A, B] = contour(Tx, Ty, C3, C3_lvl,'r');
clabel(A, B);
hold on;
title('Earth->Mars: Pork-Chop Plot');
xlabel(['X = T+ Earliest Departure Time: ', datestr(t11)]);
ylabel(['Y = T+ Earliest Arrival Time: ', datestr(t21)]);
% legend(A, 'C3 (km^2/sec^2)');



%-------------------------------------------------------------------------------

% %       [~,~,~,toh,tom,tos]...
% %          = datevec(between(TS,T1,'time'));        %Time in hours
% %       T1s = 60*(60*toh + tom) + tos;               %Time of Flight in Seconds
% %     MS = MMS * (T1s);
% %     Vo_1 = Vof(MS,e1);
% %     rval = p1/(1+e1*cosd(Vo_1));
% %       R_PQW = [rval*cos(Vo_1), rval*sin(Vo_1), 0];
% %       R_ijk = TransS.*(R_PQW);
% %         RS = [R_ijk(1,1)+R_ijk(1,2)+R_ijk(1,3),...
% %               R_ijk(2,1)+R_ijk(2,2)+R_ijk(2,3),...
% %               R_ijk(3,1)+R_ijk(3,2)+R_ijk(3,3)];          
% %    %Velocity:
% %       V_PQW = [vval1*-sin(Vo_1), vval1*(e1 + cos(Vo_1)), 0];
% %       V_ijk = TransS.*V_PQW;
% %         VS = [V_ijk(1,1)+V_ijk(1,2)+V_ijk(1,3),...
% %               V_ijk(2,1)+V_ijk(2,2)+V_ijk(2,3),...
% %               V_ijk(3,1)+V_ijk(3,2)+V_ijk(3,3)]; 

% %         [~,~,~,toh,tom,tos]...
% %            = datevec(between(TT,T2,'time'));        %Time in Hours
% %         MT = MMT * (60*(60*toh + tom) + tos);
% %         Vo_2 = Vof(MT,e2);
% %         rval = p2/(1+e2*cosd(Vo_2));
% %           R_PQW = [rval*cos(Vo_2), rval*sin(Vo_2), 0];
% %           R_ijk = TransT.*(R_PQW);
% %             RT = [R_ijk(1,1)+R_ijk(1,2)+R_ijk(1,3),...
% %                   R_ijk(2,1)+R_ijk(2,2)+R_ijk(2,3),...
% %                   R_ijk(3,1)+R_ijk(3,2)+R_ijk(3,3)]; 
% %        %Velocity:
% %           V_PQW = [vval2*-sin(Vo_2), vval2*(e2 + cos(Vo_2)), 0];
% %           V_ijk = TransT.*V_PQW;
% %             VT = [V_ijk(1,1)+V_ijk(1,2)+V_ijk(1,3),...
% %                   V_ijk(2,1)+V_ijk(2,2)+V_ijk(2,3),...
% %                   V_ijk(3,1)+V_ijk(3,2)+V_ijk(3,3)]; 