function [Trans] = TransMat(i,O,w)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    R11 =  cosd(O)*cosd(w) -sind(O)*sind(w)*cosd(i);
    R12 = -cosd(O)*sind(w) -sind(O)*cosd(w)*cosd(i);
    R13 =  sind(O)*sind(i);
    
    R21 =  sind(O)*cosd(w) +cosd(O)*sind(w)*cosd(i);
    R22 = -sind(O)*sind(w) +cosd(O)*cosd(w)*cosd(i);
    R23 = -cosd(O)*sind(i);
    
    R31 =  sind(w)*sind(i);
    R32 =  cosd(O)*sind(i);
    R33 =  cosd(i);
    
    Trans = [R11 R12 R13;...
             R21 R22 R23;...
             R31 R32 R33];


end
