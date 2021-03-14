function [V1,V2,Var,run_time,count,Lambert] = lambert(mu,Initial,Final,TOF,N,Method,IttMethod,tolerance,Lambert)
%% Solve Lambert's Problem With Several Methods
    %µ = G*M;            (mu)      [m^3/s^2]
    
    %Initial Position: (Initial) [X1m Y1m Z1m]
    
    %Final Position:    (Final)  [X2m Y2m Z2m]
    
    %Time of Flight:     (TOF)       [Sec]
    
    %Methods: 
    ...1)Universal Variable Method 
    ...2)Lagrange Method 
    ...3)Gauss Method
    ...4)Giulio Avanzini
       
    %Itt_Methods: 
    ...1)Bisection Method 
    ...2)Newton's Methond
    
    %Tolerance: 
    
    
    
%% Fixed Values:
     tic
    %Parameters: (m)
       R1 = Initial;Lanbert.R1=R1;
         r1 = norm(R1);Lanbert.r1=r1;         %Magnitude of Radius 1          (m)
       R2 = Final;Lanbert.R2=R2;
         r2 = norm(R2);Lanbert.r2=r2;         %Magnitude of Radius 2          (m)
       
       Vi = Initial;Lanbert.Vi=Vi;
       Vf = Initial;Lanbert.Vf=Vf;
       
       V1 = zeros(1,3);Lanbert.V1=V1;
       V2 = zeros(1,3);Lanbert.V2=V2;
         
       C = R2 - R1;Lanbert.C=C;                        %Chord Between R1 and R2        (m) 
         c = norm(C);Lanbert.c=c;
         
       theta = acosd(dot(R1,R2)/(r1*r2));Lanbert.theta=theta;%(??)Angle Between R1 & R2      (deg)
       ... sqrt( r1.^2 + r2.^2 - 2*r1.*r2.*cos(theta));
       s = (c+r1+r2)/2;Lanbert.s=s;           %Half Perimeter of Triangle     (m)
       am = s/2;Lanbert.am=am;                %Min Value of Semi-major axis   (m)
       ... (1/4)*(r1 + r2 + c);
       A = sqrt(r1*r2*(1+cosd(theta)));Lanbert.A=A; %Area of Triangle         (m^2)
       ... r1*r2*sin(theta)/2;               
       ... sqrt(r1*r2)/sqrt(1-cos(theta))    
       lam = (sqrt(r1*r2)*cosd(theta/2))/s;Lanbert.lam=lam;  %
 
 
%% Solver:
while abs(strcmp('finished',Method)-1)
  switch Method
  case {1,'Universal Variable','UV'}      
%% 1)Universal Variable Method
      %Equation:
         C2= @(Xi) (1 - cos(sqrt(abs(Xi))))/Xi;      Lanbert.C2=C2;
         C3= @(Xi) (sqrt(abs(Xi)) - sin(sqrt(abs(Xi))))/((sqrt(abs(Xi)))^(3/2));Lanbert.C3=C3;
         Y = @(Xi,C2,C3) r1+r2 - A*(1-Xi*C3)/sqrt(C2);Lanbert.Y=Y;
         X = @(Y,C2) sqrt(Y/C2);Lanbert.X=X;
         FXi=@(X,Y,C3) ( (C3 * X^3) + (A * sqrt(Y)) )/sqrt(mu);Lanbert.FXi=FXi;
         F = @(Xi) FXi(X(Y(Xi,C2(Xi), C3(Xi)),C2(Xi)), Y(Xi,C2(Xi), C3(Xi)),C3(Xi)) - TOF;Lanbert.F=F;
           
      %Iteration:
         switch IttMethod
         case {1,'Bisection','bisection','bi','Bi','B','b'}
             [Xi,count] = bisection(F,am,-20,tolerance);
         case {2,'Newton','newton','newit','Newit','N','n'}
             [Xi,count] = newit(F,am,tolerance);  
         otherwise
             warning(...
             'No Method was picked. Program will default to bisection Method');
             [Xi,count] = bisection(F,am,-20,tolerance);
         end
         Lanbert.Xi=Xi;Lanbert.count=count;
         
      %Finding Velocity Vector(s):    
         c2= C2(Xi);Lanbert.c2=c2;
         c3= C3(Xi);Lanbert.c3=c3;
         y = Y(Xi,c2,c3);Lanbert.y=y;
         f = 1 - y/r1;Lanbert.f=f;
         g = A*sqrt(y/mu);Lanbert.g=g;
         V1= (R2 -f*R1)/g;Lanbert.V1=V1;
         gp= 1 - y/r2;Lanbert.gp=gp;
         fp= (sqrt(mu*y)*(-r1 -r2 +y))/(r1*r2*A);Lanbert.fp=fp;
         V2 = fp*R1 + gp*V1;Lanbert.V2=V2;
         
       Var = Xi;    Lambert.Var=Var;
       Lambert.Variable = 'Xi';
       
       Method = 'finished';
            
  case {2,'Lagrange','L'}
%% 2)Lagrange Method
      %Equation:
         fy = @(x) sqrt(1 - (lam^2) * (1-x^2));Lanbert.fy=fy;
         fz = @(x,y) y - lam*x;Lanbert.fz=fz;
         fS1= @(x,y,z) (1 - lam - x*z)/2;Lanbert.fS1=fS1;
         fQ = @(x,S1) (4/3) * hypergeom([3,1],5/2,S1);Lanbert.fQ=fQ;
         fTOF = @(x,z,Q) ((z^3)*Q + 4*lam*z)/sqrt(mu/(am^3));Lanbert.fTOF=fTOF;
         f = @(x) fTOF(x,fz(x,fy(x)),fQ(x,fS1(x,fy(x),fz(x,fy(x)))));Lanbert.f=f;
          
         alpha = @(a)...
             (2*asin(sqrt(s/2*a)));
          ...(2*asin(sqrt((a^3)/mu)));
              Lanbert.alpha=alpha;
         beta = @(a)...
             (2*asin(sqrt((s-c)/(2*a))));
             Lanbert.beta=beta;
            
         FA = @(a)...
             sqrt((a^3)/mu)*((alpha(a)-sin(alpha(a)))-(beta(a)-sin(beta(a))))-TOF;
             Lanbert.FA=FA;
         
         
         
      %Iteration:
         switch IttMethod
         case {1,'Bisection','bisection','bi','Bi','B','b'}
             [Var,count] = bisection(FA,am,-20,tolerance);
         case {2,'Newton','newton','newit','Newit','N','n'}
             [Var,count] = newit(FA,am,tolerance);  
         otherwise
             warning(...
             'No Method was picked. Program will default to bisection Method\n');
             [Var,count] = bisection(FA,am,-20,tolerance);
         end
         a = Var;
         Lanbert.a=a;
         Lanbert.Var=Var;
        Lambert.Variable = 'a';
             
      %Finding Velocity Vector(s):  
         phi = (alpha(a) + beta(a))/2;Lanbert.phi=phi;
         psi = (alpha(a) + beta(a))/2;Lanbert.psi=psi;
         x = sqrt(1 - (am/a));Lanbert.x=x;
         UVR1 = R1/(r1);Lanbert.UVR1=UVR1;
         R1_ih= [UVR1(1) 0 0];Lanbert.R1_ih=R1_ih;
         H_ih = cross(R1,R2)/norm(cross(R1,R2));Lanbert.H_ih=H_ih;
         y = fy(x);Lanbert.y=y;
         z = fz(x,y);Lanbert.z=z;
         fV1 = @(x) (1/z)*(sqrt(mu/am))*((2*lam*(am/r1)-(lam + x*z))*...
             R1_ih+sqrt(r2/r1)*(1/2)*theta*pi/180*cross(H_ih,R1_ih));Lanbert.fV1=fV1;
         V1 = fV1(x);Lanbert.V1=V1;
         
       Method = 'finished';
             
  case {3,'Gauss','G'}
%% 3)Gauss Method
      %Equation:
         tm = 1;Lanbert.tm=tm;      % +1 or -1
         costheta = dot(R1,R2)/(r1*r2);Lanbert.costheta=costheta;
         sintheta = tm*sqrt(1-(cosd(theta)).^2);Lanbert.sintheta=sintheta;
         L = ((r1 + r2) / (4* sqrt(r1*r2) * cosd(theta/2))) - (1/2);Lanbert.L=L;
         m = (mu*(TOF.^2)) / ((2* sqrt(r1*r2) * cosd(theta/2))^3);Lanbert.m=m;
          
      %Iteration:
         switch IttMethod
         case {1,'Bisection','bisection','bi','Bi','B','b'}
              [y,count] = bisection(g,-20,0,tolerance);
         case {2,'Newton','newton','newit','Newit','N','n'}    %Initial Guess:(x0,k)
%             [y,count] = newit(g,0,tolerance);
              runs = 10;
              err = zeros(1,runs);
              y = zeros(1,runs);
              y(1) = 10;
              for i = 1:runs
                  x1 = (m/(y(i)^2))- L;        %First equation of Gauss
                 %X2:
                    x2 = 1;
                    for j = 1:100
                        x2 = x2 + x2*(((4+2*j)*x1)/(3+2*j));
                    end
                  x2 = (4/3) * x2;
                  y(i+1)= 1 + x2 * ( L + x1 );
                  err(i)=abs((y(i+1)-y(i))/y(i));
                  if err(i)<tolerance
                     break
                  end
              end
              count = i;
         end
         
         Lanbert.y=y;
         Var = y(i);Lanbert.Var=Var;
            
         Lanbert.err=err;Lanbert.runs=runs;
         Lanbert.x1=x1;Lanbert.x2=x2;
         Lanbert.count=count;
                
      %Finding Velocity Vectors:  
         f = 1 - y(i)/r1;     Lanbert.f=f;
         g = A*sqrt(y(i)/mu); Lanbert.g=g;
         V1 = (R2 -f*R1)/g;   Lanbert.V1=V1;
          
         gp = 1 - y(i)/r2;    Lanbert.gp=gp;
         fp = (sqrt(mu*y(i))*(-r1 -r2 +y(i)))/(r1*r2*A);Lanbert.fp=fp;
         V2= fp*R1 + gp*V1;   Lanbert.V2=V2;
         Lambert.Variable = 'y';
         
       Method = 'finished';
  case {4,'Giulio Avanzini','GA'}      
%% 4)Giulio Avanzini
      %Assign Problem Geometry & Transfer Time:	(?,??,T12s) 
        %?  = rho
         rho =  r2/r1; Lanbert.rho=rho;                        % ratio of radii
        %T12s=TOF
         n_ref = sqrt(mu/(r1^3));
         t_ref = 1/n_ref;
         TOF_nds = n_ref * TOF;
        %?? = theta
               
         e_max = -1/cosd(theta/2);
         e_F = (r1 - r2)/c;
         e_P = sqrt(1-e_F^2);
         e_H = sqrt(e_max^2 -e_F^2);
              
      %Is ?? > 0? & Derive Limit For e_t:        (e_H,e_P)
         upper_lim = e_P;
         lower_lim = -1e20;
         if theta > 0              %Case 3 & 4
            lower_lim = -1 * e_H;
            eTf = @(x) e_P * e_H * ((x-1)/(e_P + e_H*x));
         else          %Case 1 & 2
            eTf = @(x) e_P * (1 - exp(-x/e_P));
         end
          
         g = @(x) log(hypergeom(theta*pi/180,rho,eTf(x)));
          
      %Itteration: 
         switch IttMethod
         case {1,'Bisection','bisection','bi','Bi','B','b'}
              [e_T,count] = bisection(g,upper_lim,lower_lim,tolerance);
           
         case {2,'Newton','newton','newit','Newit','N','n'}    %Initial Guess:                            (x0,k) 
%             [e_T,count] = newit(g,0,tolerance);
              runs = 20+1;
              err = zeros(1,runs);
              y = zeros(1,runs);
              y_s = log(TOF_nds);
              x = zeros(1,runs);
              x(1) = 0;   
              k = 0+1;
              Dif = diff(sym(g));
              DIF= matlabFunction(Dif);
              while 1 
      %Inverse Transform:
                    e_T = eTf(x(k));
      %Kepler's Time Equation:
                    TOF_ndg = hypergeom(rho,theta*pi/180,e_T);
                    y(k) = log(TOF_ndg);
      %Logic Test 1:
                    err(k)= abs( y(k) - y_s );
                    if err(k)<tolerance
                       break
                    end
      %Logic Test 2:
                    if lower_lim > e_T || e_T > upper_lim
                       warning('Out of Bounds');
                       break
                    end
      %Logic Test 3:
                    if k > runs
                       warning('No Convergence');
                       break
                    end
      %Incremented Unknown::
                    dx = 1;
                    for j = 1:100
                        dx = dx + dx*(( (4+2*j)*x(k)) / (3+2*j) );
                    end
                    dx = (4/3) * dx;
                    xp = x(1) + dx;
      %Inverse Transform:
                    e_Tp = eTf(x(k));
      %Kepler's Time Equation:
                    TOF_ndp = hypergeom(rho,theta*pi/180,e_Tp);
                    y(k+1) = log(TOF_ndp);
      %N-R Iteration:
                    dydx1 = (y(k+1) - y(k))/dx;
                    dydx2 = real(DIF(x(k)));
                    x(k+1) = x(k) - (y_s - g(x(k)))/dydx1;
                    k = k + 1;
              end
              count = k;
         end      
               
      %Finding Velocity Vectors:
         wc = acosd(dot(R1,C)/(r1*c));      Lanbert.wc=wc;
         a_F = (r1 + r2)/2;                 Lanbert.a_F=a_F;
         p_F = a_F * (1 - e_F);             Lanbert.p_F=p_F;
         e = sqrt(e_F^2 +e_T^2)-.2;         Lanbert.e=e;
         p = p_F - e_T *((r1*r2)/c)* sind(theta);Lanbert.p=p;
         a = p/(1-e^2);                     Lanbert.a=a;
         EV= [(e_F*cosd(wc) - e_T*sind(wc)),(e_F*sind(wc) + e_T*cosd(wc)),0]; Lanbert.EV=EV;
         w = atand((e_F*sind(wc) + e_T*cosd(wc))/(e_F*cosd(wc) + e_T*sind(wc)));Lanbert.w=w;
         Vo = -w; Lanbert.Vo=Vo;
         Vof = theta - w; Lanbert.Vof=Vof;
         i = 0; Lanbert.i=i;
         Lanbert.Omega = -60;
%        Va = @(Vb,Vc)  -(R(1)*mu*e - (R(1)*r1*(Vb^2)) - (R(1)*r1*(Vc^2)) )/...
%                       (r*(Vb*R(2) + R(3)*Vc));
%        Vb = @(Va,Vc)  -(R(2)*mu*e - (R(2)*r1*(Vc^2)) - (R(2)*r1*(Va^2)) )/...
%                       (r*(Vc*R(3) + R(1)*Va));
%        Vc = @(Vb,Va)  -(R(3)*mu*e - (R(3)*r1*(Va^2)) - (R(3)*r1*(Vb^2)) )/...
%                       (r*(Va*R(1) + R(2)*Vb));
           
         V1 = [r1*cosd(Vo) r1*sin(Vo) 0];   Lanbert.V1=V1;
         V2 = [r2*cosd(Vof) r2*sin(Vof) 0]; Lanbert.V2=V2;
         Var = e_T;              Lanbert.Var=Var;
%          Var = [a e 0 0 w Vo];              DATA.Var=Var;
         Lambert.Variable = 'e_T';
           
        Method = 'finished';
        
  case {5}      
%% 5)
        E6i(1:3) = R1;  E6i(4:6) = Vi;
        E6f(1:3) = R2;  E6f(4:6) = Vf;
        [V1,V2]= newlambert(mu,E6i,E6f,TOF,N);
        Var = 1;
        count = 0;
       Method = 'finished';
  otherwise
%%  X)No Method Picked
      warning(...
       'No Method was picked. Probram will default to Universal Variable Function\n');
      Method = 1;
  end    
end    
  run_time = toc;Lanbert.run_time=run_time;             %Run Time Of Method   (Sec)
  Lambert.check = Var;
end
