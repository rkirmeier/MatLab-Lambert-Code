function [x,counter] = bisection(Fun,a,b,t)
%Bisection Iterative Method:
% fx = equation with variable x
    digits (16);
    format long;
     counter = 0;
     Max_Runs = 1000;
      x = (a+b)/2; 
%       Fx= Fun(x);
      Fa= Fun(a);
      Fb= Fun(b);

       if Fa * Fb > 0      % Display error and finish if signs are not different
          fprintf(1,'Have not found a change in sign. Will not continue...\n');
          x = 0;%'ERROR';
          return
       end
       while(abs(b-a)>t)
          x=(a+b)/2;
          Fx=Fun(x);
          if((Fx==0) || (abs(Fx)<t)) %If F(x)=0, x is a root
            a=x; b=x;    %If |F(x)|< tollerance, x is a root
          else
              if(Fb*Fx>0)
                 b=x; Fb=Fx;
              end
              if (Fa*Fx>0)
                 a=x; Fa=Fx;
              end
          end
          counter = counter + 1;
          if counter > Max_Runs
             error('Does not converged after %5f iterations\n',counter); 
             return
          end
      end
end

