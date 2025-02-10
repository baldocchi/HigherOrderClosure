function [can] = CANLEN(can)



    
% Compute canopy length scale.  The constant C0 is used in these computations

      
      c0d=can.C0/can.drag;
	  
       can.l(1,1)=0.;

	   nlay=can.sze1+1;

	   % test if dl/dz <= k, (0.40)


% if expression
%     statements
% elseif expression
%     statements
% else
%     statements
% end


      for I=2:nlay
	  
        if(can.pad(I,1) < 0.0001)
		
           can.l(I,1)=can.l(I-1,1)+can.grad; 
			can.ld(I,1)=can.l(I,1)*can.gamma;
       
        elseif (can.pad(I,1) >= 0.0001)
		
		   
         can.l(I,1)=c0d/can.pad(I,1);       % need to divide by LAI density to create length scale

         if(can.l(I,1)-can.l(I-1,1) > can.grad)
			 can.l(I,1)=can.l(I-1,1)+can.grad ;
         end
         
		if(can.l(I,1)-can.l(I-1,1) < -can.grad)
			can.l(I,1)=can.l(I-1,1)-can.grad;
        end
		   
        can.ld(I,1)=can.l(I,1)*can.gamma;
		
        end 
	  
     

end