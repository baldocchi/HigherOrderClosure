function [can] = LGEN(can)

      
      nm1=can.sze3-1;
      
	  for I=nm1/3+2:can.sze3
	  
         can.l(I,1)=can.l(I-1,1)+.02;
         can.ld(I,1)=can.l(I,1)*can.gamma;    % dissipation length scale
	 end
end