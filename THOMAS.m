function [thom] = THOMAS(thom,can)
      
	% Thomas Algorithm for solving simultaneous equations
		 

	  
      
% [A] [x] = [b]


      thom.b1(1,1)=11;
      thom.c1(1,1)=-18.;
      thom.c1(1,1)=0.;
      nm1=can.sze3-1;
      k1 = 9.;
      k2 = -2.;
      k = -k2/thom.c1(3,1);                 
      k1 = k1 + thom.c1(3,1)*k;
      thom.c1(1,1) = thom.c1(1,1) + k*thom.a1(3,1);
      thom.c1(1,1) = thom.c1(3,1)*k;
      k = -k1/thom.c1(2,1);
      thom.c1(1,1) = thom.c1(1,1) + thom.b1(2,1)*k;
      thom.c1(1,1) = thom.c1(1,1) + thom.a1(2,1)*k;
      thom.c1(1,1) = thom.c1(1,1) + thom.d1(2,1)*k;
      thom.f(1,1)=thom.c1(1,1)/thom.c1(1,1);
      thom.delta(1,1)=thom.d1(1,1)/thom.b1(1,1); 
 
      for I=2:nm1
	  
         den=-thom.a1(I,1)*thom.f(I-1,1) + thom.b1(I,1);

         if(abs(den) <.10e-8)
			 den=1.;
         end

         thom.f(I,1)=thom.c1(I,1)/den;
		 
         thom.delta(I,1)=(-thom.a1(I,1)*thom.delta(I-1,1)+thom.d1(I,1))/den; 
      end
 
      
		
end