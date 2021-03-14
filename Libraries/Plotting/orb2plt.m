function [body] = orb2plt(Dimension,Size,Color,Orbit,Position,Velocity)
%UNTITLED3 Summary of this function goes here
%% *History*
%  When       Who    What
%  ---------- ------ --------------------------------------------------
%  2019/07/11 mnoah  original code
%  2020/21/10 kirmeier  the butcherer of original code

%% Plot:
switch Dimension
    case {1,2}
      %Detailed explanation goes here  
            p = Position(1);  q = Position(2);  % [m] Position
            pv = Orbit(1);    qv = Orbit(2);    % [m] Orbit
% 
%             p1 = pq1(1);    q1 = pq1(2);% [m] Velovity Arrow
%             xi = xy{1};     yi = xy{2}; % [m] Velovity Arrow (inside)
%             xr = xy{3};     yr = xy{4}; % [m] Velovity Arrow (outside)
       %Plot Orbit:  
         plot3(pv,qv,0,'-b');
   
       %Plot Position: 
%          plotbody(Dimension,Position,Size,0.01,Color);%'.r'
         %plot(p,q,'.r','MarkerSize',25);
   
%        %Plot Velocity Arrow: 
%          plot([p p1],[q q1],'r');
%            plot(xi,yi,'r');
%            plot(xr,yr,'r');
           
    case 3   
        
       %Plot Orbit:  
          plot3(Orbit(:,1),Orbit(:,2),Orbit(:,3));
       %Plot Position: 
          body = plotbody(Dimension,Position,Size,36,Color);%'.r'

         %plot(p,q,'.r','MarkerSize',25);
   
%        %Plot Velocity Arrow: 
%          plot([p p1],[q q1],'r');
%            plot(xi,yi,'r');
%            plot(xr,yr,'r');
    otherwise
        fprintf("Pleace Pick A dimension Between 1 & 3 and try again")
        return
        
end
     
% % Axies
%    plot([-1.25*rp 1.25*rp],[0 0],'k');
%    plot([0 0],[-1.25*ra 1.25*ra],'k');
%    text(1.25*rp, 0.1*rPlanet, 'p');
%    plot(1.25*rp+[0 -0.1*rPlanet*cosd(30)],[0 0.1*rPlanet*sind(30)],'k');
%    plot(1.25*rp+[0 -0.1*rPlanet*cosd(30)],[0 -0.1*rPlanet*sind(30)],'k');
%    text(0.1*rPlanet, 1.25*ra, 'q');
%    plot([0 -0.1*rPlanet*sind(30)],1.25*ra+[0 -0.1*rPlanet*cosd(30)],'k');
%    plot([0 +0.1*rPlanet*sind(30)],1.25*ra+[0 -0.1*rPlanet*cosd(30)],'k');


end

