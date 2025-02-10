function [can] = ZPD(turb,can)
 
	% Zero Plane Displacement

    % Thom's center of pressure method is used, such that d is associated
    % with the height where 1/2 of the momentum transferred to the surface is
    % absorbed.

    
     
      sfmdg=0.0; 
      zfmdg=0.;
      fz1=0.;
      f1=0.0;
      j=2 ;
      nm1=can.sze3-1;
	  grd=1./(nm1/3) ;
	  zz=grd;
      nlyr=nm1/3+1;  % number of layers
      
	  for I=1:nm1/6
	     f2=-turb.duw(j,1); 
         f3=-turb.duw(j+1,1); 
         fz2=f2*zz;
         zz = zz+ grd;
         fz3=f3*zz;
         sfmdg = sfmdg+ grd/3.*(f1+4.*f2+f3);
         zfmdg =  zfmdg + grd/3.*(fz1+4.*fz2+fz3);
         f1=f3;
         fz1=fz3;
         zz = zz+ grd;
         j = j+ 2;
       end

	  % zero plane displacement

      can.d=zfmdg/sfmdg*(turb.uw(nlyr+1,1)-turb.uw(1))/turb.uw(nlyr+1,1);
end

      
  
 
