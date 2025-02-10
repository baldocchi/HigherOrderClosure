function [can]=GENLD(L,AL,NG)
     for I=1:NG
         can.LD(I)=L(I)*AL;
     end
end
     