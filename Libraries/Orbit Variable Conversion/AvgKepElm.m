function [Element] = AvgKepElm(Initial,Final)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

 Element.R = Initial.R;
 Element.V = Initial.V;
 Element.a = (Initial.a + Final.a)/2;
 Element.e = (Initial.e + Final.e)/2;
 Element.i = (Initial.i + Final.i)/2;
 Element.Omega = (Initial.Omega + Final.Omega)/2;
 Element.w = (Initial.w + Final.w)/2;
 Element.Vo_1 =Initial.Vo;
 Element.Vo_2 =Final.Vo;
 Element.Extra =' ';
 Element.H = (Initial.H + Final.H)/2;
 Element.N = (Initial.N + Final.N)/2;
 Element.p = (Initial.p + Final.p)/2;
 Element.EV =(Initial.EV + Final.EV)/2;
 Element.II =(Initial.II + Final.II)/2;
 Element.ra = (Initial.ra + Final.ra)/2;
 Element.rp = (Initial.rp + Final.rp)/2;
 Element.TP = (Initial.TP + Final.TP)/2;
 Element.MM = (Initial.MM + Final.MM)/2;
 Element.Ei = Initial.E;
 Element.Mi = Initial.M;
 Element.Lo1 =Initial.Lo;
 Element.uo1 =Initial.uo;
 Element.E2 = Final.E;
 Element.M2 = Final.M;
 Element.Lo2 =Final.Lo;
 Element.uo2 =Final.uo;
 Element.Vo =Initial.Vo;
end

