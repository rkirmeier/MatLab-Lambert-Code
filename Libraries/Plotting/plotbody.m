function [body] = plotbody(Dimension,Position,radius,tollerance,Color)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
   Evals = 0:tollerance:360.0; % [deg] values of the eccentric anomaly around orbit

   switch Dimension
       case 1
         body = plot(Position(1),Position(2),Color,'MarkerSize',radius);    %plot(,Position(2),'o'); 
       case 2
         body = fill(Position(1)+radius*cosd(Evals),Position(2)+radius*sind(Evals),Color);
%           body = fill(Position(1) + RX,Position(2) + RY, Color);
       case 3
         [X,Y,Z] = sphere(360/tollerance);
            RX = X*radius;
            RY = Y*radius;
            RZ = Z*radius;
           body = mesh(Position(1)+RX, Position(2)+RY, Position(3)+RZ); %,'FaceColor',Color);
           body.FaceColor = Color;
           body.EdgeColor = Color;
       otherwise
           body = 'Error';
           return
   end
end

