function [root,i] = newit(Fun,g,t)
% Program Code of Newton-Raphson Method in MATLAB 
% Fun = equation with variable x 
% Fun=input('Enter the function in the form of variable x:','s')
% x(1)=input('Enter Initial Guess:')
% t=input('Enter allowed Error:');
     runs = 100000;
     err = zeros(1,runs);
     x = zeros(1,runs);
      x(1) = g;
%      f=inline(Fun);
%      DIF=diff(sym(Fun));
     Dif = diff(sym(Fun));
      DIF= matlabFunction(Dif);
     for i=1:runs
        x(i+1)=x(i)-((Fun(x(i))/DIF(x(i))));
        err(i)=abs((x(i+1)-x(i))/x(i));
        if err(i)<t
            break
        end
     end
     root=x(i);
end




