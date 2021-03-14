function [IJK] = matsolv(PQW,TRANS)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

      %R_IJK = R_Trans*R_PQW;
          ijk = TRANS.*PQW;
            IJK(1) = ijk(1,1)+ijk(1,2)+ijk(1,3);
            IJK(2) = ijk(2,1)+ijk(2,2)+ijk(2,3);
            IJK(3) = ijk(3,1)+ijk(3,2)+ijk(3,3);  
end

