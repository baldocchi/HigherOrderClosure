% HigherOrderClosure_v2025.m : 
% 
% Matlab version of the C version of the code I made of the Fortran version 
% of Meyers and Paw U Higher Order Closure model
% 
% 
%  Meyers TP + Paw U KT. 1987. Ag.Forest Met.
% 
% Meyers , T. P., and K. T. Paw U (1987), 
% Modeling the plant canopy micrometeorology with higher order 
% closure principles, Agricultural and Forest Meteorology, 41, 143-163.

%  
%  Dennis Baldocchi
%  Ecosystem Science Division
%  Dept of Environmental Science, Policy and Management
%  University of California, Berkeley
%  baldocchi@nature.berkeley.edu


% Feb 5, 2025

% Tilden gave us an updated version of Fortran in 2025

%  Dec. 2003 (From Tilden's 1987 Fortran version)
% 
%  Using code to compute turbulence statistics in sparse savanna canopy
% 
%  Nine equations are solved for
% 
%  <u>, <w'u'>,<u'u'>, <v'v'>, <w'w'>,<w'u'u'>, <w'v'v'>, <w'w'w'>, <w'w'u'>
% 
%  Fourth order terms are solved with the quasi-Gaussian approximation
% 
%  <w'w'w'w'>=3 <w'w'>^2
% 
%  <w'w'w'u'> = 3 <w'w'><w'u'>
% 
%  <w'u'u'u'> = 3 <u'u'> <w'u'>
% 
%  <w'w'u'u'> = <w'w'><u'u'> + 2 <w'u'><w'u'>
% 
%  <u'u'u'u'> = 3 <u'u'><u'u'>
% 
%  Computations are made for a normalized canopy with x layers

%  The grid domain of the canopy extends to three times canopy height times layers
% 
%  The turbulence profiles are also normalized in terms of friction velocity, u*
% 
%  the model assumes stationarity



%  Debugging

% Feb 2025 

% increased Cd to 0.2 and get a 2nd wind max!!

% refining the tests for convergence and it is looking better

% the u profile seems to be approaching 0 in the canopy. need to fix a few
% more things

% Tilden's Cionco code U function was extinquishing u in the canopy too,
% too fast. I needed to normalize the layers by the number of canopy layers

% was still having problems and getting nonsense U profiles.  then notice
% code did not mimic the numeric for solving the equation. made change and
% now profiles seem to better follow LAI profile

% Feb 5 I went and looked at the Adams Bashford integration and implemented
% that for u profiles. seems to work better, and using trapezoid for layer
% 2

%  Jan 2025

% with matlab I can plot turbulence profiles and with each iteration.  It
% seems that the initial conditions of w2,u2 and v2 follow a linear
% function of uw.  When I plot the iterations these values are not
% changing. They dont seem to reflect the LAI I am giving them, like this
% case with Harvard forest and a top heavy LAI.

% just set the wt at 0.25 and I get some interesting profiles and lots of
% iterations

% At present the only profile that looks funky is wind velocity in the
% canopy

%  so dive into that code


% Plots are showing the upper nodes being noisy and xx(sze3) not zero or
% one..


% few more debugs.. some noise at top of domain for xxx and kink at the top
% of canopy for wu
% not seeing secondary wind max in the canopy.

% 
%  June 2010...found some old notes stating that James found the code to be sensitive to number of layers
%  just ran 100 layer model and found u* normalized to 1 finally!! Also d/h is 0.79 and zo/h is 0.09
%  reasonable numbers??..though profiles of the variances are pretty linear..not much action from the canopy
%  I also got lower d/h (0.63) when I reduced the lai profile...so this is a good sign
%  tried to fix iterations..set old values to 0.. It was equal to U, but then the interations were on qq
%  last run had 79 iterations..A chance to move from the initial conditions.  Early runs were stopping
%  after less than 10 iterations.  I also reset iterations to test only inside the canopy..The above canopy
%  profile was invariant and was biasing the 'fit' and causing the model to 'converge' too soon. Now 
%  I am up to several hundred iterations


% % found some errors in layers of canlen and passing variables in and out of THOMAS
% 
% % check ZPD more
% 
% % Dec 14, 2003  The program seems to be working. Computations of profiles produced
% 
% % Recasted the interation do loop into a while statement
% 
% % Jan. 5, 2004
% % Fixed lai profile to scale non-dimensionally
% 
% % Jan 7, 2004
% 
% % debug runs are showing a strong sensitivity of turbulence profiles to LAI
% % profiles. Getting kinks at interface of ww and uu profiles if lai is not
% % skewed and peaked close to the canopy-atmopshere interface. Solutions also seem 
% % to blow up
% 
% % Changed convergence criterion. Expanded itermax to 1500 and am converging with
% % the integration of the wind profile (prior and current) to 10-3. Also made
% % the weighting function weights more strict.
% 
% % iterations number around 550 now and smoother profile
% 
% 
% %Jan 9, 2004
% % There seems to be some inconsistencies. Many calculations are done with u/u* and w'u'/u*^2
% % But the variance and triple products seem to be computed with normal profiles, then they
% % and u/u* are normalized by u* again. I was worried because u/u* was near 12, which is to large.
% % Now it is near 6 at the top of the canopy.
% 
% 
% % Test of the 'normalized' output shows that u is already normalized by u*
% % at the end of the run, but we need to normalize the other turbulence statistics
% % in order to  match the initial normalized boundary conditions...search continues
% 
% 
% % as iterations continue wu drops from -1 to about -0.8 or so. Part of issue seems
% % to be feedback on dwu/dz ~Cd a u u
% 
% % If I set test for sumu at 0.01 (out of 500.00 or one out of 10,000) I then get
% % convergence after 7 iterations or so and the 'real' ustar is near 1.000
% % which makes sense, since we are working with in terms of a normalized u*
% % system anyhow.  There are no initial conditions for u or u*
% 
% % plus sigw/U* is about 1.25 and the other coefficients make sense
% 
% % Jan 11 model produces too much ww and uu in upper canopy. Try alternative
% % length scale L~ du/dz/uh
% 
% 
% % March 19, 2005. Computing vertical profile of v to examine veering of wind directino
% % profile within the canopy
% 
% % April 21, 2005 I was getting different solutions and number of iterations whether
% % I tested for u, qq or a higher order moment like wwu.  In particular I was getting
% % a few iterations and u* of near 1 with u and 200+ iterations and a u* 0f 0.5
% % with qq.
% 
% % Then I reset the filter coefficients and kept them strict, eg wt =0.01
% % Now I get u* near 1 for all three classes (u, qq, wuu) of convergence criteria
% % but wwu yields about 200 iterations and qq about 21.
% 
% % use paired t test to iterate
% 
% % Oct 4, 2006
% 

clear all;
close all;

tic;
can.sze=100;
can.sze1= 101;					% the number of canopy layers is 101 
can.sze3= 301;				    % number of grid points, 3 times canopy layers plus one

can.C0 = 0.060;    % constant used to define turbulence length scale
C2 = 3.62;	   % higher order closure coefficients for second order terms
C3 = 6.92; 
C8 = 10.5;     % higher order closure coefficients for third order terms

ALPH0= 0.15;  % constants in the trace free component of pressure term <u'dp'>/dx
ALPH1= 0.43;	 
ALPH2 = 0.23;	

% In Meyers and Paw U alph0=0.30


dthom = 0.60;  % zero plane displacement, Thom center of pressure method, height of 0.5<w'u'>; as a function of canopy height
GAMMA = 18;    % constant applied to dissipation length scale, e.g. epsilon = q^3/ GAMMA l
can.z0=0.1;    % roughness length as a function of canopy height




% Declare Structures

% in Matlab assign structures as zeros

% structure declaring turbulence variables and arrays

% first order terms

turb.u=zeros(can.sze3,1); % mean horiziontal wind velocity, m/s
turb.us=zeros(can.sze3,1);% mean wind speed, us=u(1+<v'v'>/2<u><u>)
turb.v=zeros(can.sze3,1);% mean v velocity
turb.du=zeros(can.sze3,1); % wind shear, du/dz
turb.ws=zeros(can.sze3,1); % wind shear, du/dz
turb.w=zeros(can.sze3,1);   % <w>
turb.uspeed=zeros(can.sze3,1); % sqrt(u2 + <u'u'> + v2 + <v'v'> + w2 +<w'w'>)
turb.mw=zeros(can.sze3,1); % mean horiziontal wind velocity, m/s


% second order terms

turb.u2=zeros(can.sze3,1);   % longitudinal velocity variance
turb.v2=zeros(can.sze3,1);   % lateral velocity variance
turb.w2=zeros(can.sze3,1);   % vertical velocity variance
turb.uw=zeros(can.sze3,1);  % Reynolds shear stress
turb.ustar=zeros(can.sze3,1) ;     % friction velocity
turb.duw=zeros(can.sze3,1); % d<w'u'>/dz

% third order terms	 

turb.dwwu=zeros(can.sze3,1); % d<w'w'u'>/dz
turb.wwu=zeros(can.sze3,1);  % <w'w'u'>
turb.w3=zeros(can.sze3,1);   % <w'w'w'>
turb.dw2=zeros(can.sze3,1);  % d<w'w'>/dz
turb.wuu=zeros(can.sze3,1);  % <w'u'u'>
turb.wvv=zeros(can.sze3,1);  % <w'v'v'>
turb.wuu1=zeros(can.sze3,1); % <w'u'u'> term without du/dz 
turb.dwuu1=zeros(can.sze3,1); % d<w'u'u'1>/dz
turb.du2=zeros(can.sze3,1);	 % d<u'u'>/dz
turb.dv2=zeros(can.sze3,1);   % d<v'v'>/dz
turb.wqq=zeros(can.sze3,1);	 % <w'q'q'>
turb.dwqq=zeros(can.sze3,1);  % d<w'q'q'>/dz
turb.ust=zeros(can.sze3,1);   % friction velocity profile, ustar
turb.wvv1=zeros(can.sze3,1);  % <w'v'v'> term without du/dz 
turb.dwvv1=zeros(can.sze3,1); % d<w'v'v'1>/dz
turb.w31=zeros(can.sze3,1);   % <w'w'w'> term without du/dz 
turb.dw31=zeros(can.sze3,1); % d<w'w'w'1>/dz


% TKE

turb.qq=zeros(can.sze3,1);	 %  tke, qq=[u'u' + v'v' + w'w']
turb.q=zeros(can.sze3,1);     %  sqrt qq

% time scales

turb.tau=zeros(can.sze3,1);   % turbulence time scale, q^2/epsilon

% iteration test

turb.old=zeros(can.sze3,1);
	  

%  structure with information on canopy architecture

can.pad=zeros(can.sze3,1); % profile of plant area density, m2 m-3
can.l=zeros(can.sze3,1);   % length scale
can.ld=zeros(can.sze3,1);  % dissipation length scale
can.lai_z=zeros(can.sze3,1);

	 

     % otherstuff
	 	 
fact.x1=zeros(can.sze3,1);
fact.x2=zeros(can.sze3,1);
fact.x3=zeros(can.sze3,1);
fact.x4=zeros(can.sze3,1);
fact.x5=zeros(can.sze3,1);
fact.x6=zeros(can.sze3,1);
fact.x7=zeros(can.sze3,1); 
	 
fact.von_kar=zeros(can.sze3,1);  % von Karman's constant divided by u*
fact.denom=zeros(can.sze3,1);
	 
fact.x=zeros(can.sze3,1);
fact.dx=zeros(can.sze3,1);
	


% matrices passed in and out of thomas algorithm

     thom.a1=zeros(can.sze3,1);
	 thom.b1=zeros(can.sze3,1);
	 thom.c1=zeros(can.sze3,1);
	 thom.d1=zeros(can.sze3,1);
	 thom.delta=zeros(can.sze3,1);
	 thom.f=zeros(can.sze3,1);



     % start of Main

	  nm1=can.sze3-1;
      nm2=can.sze3-2;
      nm3=can.sze3-3;
      nlyr=nm1/3 + 1;    % number of canopy layers
      grid=1./double(nm1/3) ;   % grid is the layer of canopy heights
      gridsq=grid*grid ;
	  gridp2=grid+grid;
      
      


	  % input Canopy Height and Drag 
	 
      % Generate plant area density profile 

		% Will need to compute normalized LAI profile

	   [can]= LAI_PROFILE(can);
      % [can]= LAI_PROFILE_HARVARD(can);

       % [can]= LAI_PROFILE_Tonzi(can);
	        
        % Canopy info is in LAI_PROFILE subroutines

		
		can.drag=0.2;		   % canopy drag coefficient, 0.1
		can.grad=.4*grid;      % k dz, to test if dl/dz <= k (0.40)
       can.zg=.005/can.ht;    % non-dimensional roughness parameter at the soil
		%can.zg=0.01;  % soil roughness
       can.gamma=GAMMA;
        can.d=dthom;


		

 
%  Compute within and above canopy length scales


     [can]= CANLEN(can); 

     [can]= LGEN(can); 

     
 % initialize normalized log wind profile above the canopy 

% u(z)/u* = 1/k log ((z-d)/z0)

% (du/dz)/u* = 1/k(z-d)

% note 2.5 = 1/0.4 = 1/k


       % above the canopy

	  % I think u should be normalized later.

      for I=(nlyr+1):can.sze3
	  
       turb.u(I,1)=2.5*log((double(I-1)*grid-can.d)/can.z0);
       turb.du(I,1)=2.5/(double(I-1)*grid-can.d);


       % computed on a grid basis, for a normalzied canopy height

       % u(z) =u*/k ln((z-d)/z0) d=0.6h, z0=0.1h

      
        % In TDM code

        turb.mw(I,1)=2.5*log((double(I-1)*grid-can.d)/can.z0);
        turb.ws(I,1)=2.5/(double(I-1)*grid-can.d);  % wind shear gradient..check if I need to divide by dz

        end


        UH=turb.u(nlyr+1);
	
% initialize wind profile and wind shear inside the canopy using
% exponential model of Cionco u(z) = uH exp(a (z/H -1))

% 3.5=uh/u* at the top of the canopy


% i should compute Uh first from log profile instead of assume it is 3.5

      for I=2:nlyr
	  
         % turb.u(I,1)=UH*exp(2.*(double(I-1)*grid - 1.));  % estimate of wind speed in canopy using Cionco Model 
		 % turb.old(I,1)=0;
         % 
         turb.u(I,1)=UH*exp(2.*(double(I-1)/can.sze1 - 1.));  % estimate of wind speed in canopy using Cionco Model 
		 turb.old(I,1)=0;

         % in TDM code refer to mean wind and wind shear. I use u and du
         % turb.mw(I,1)=UH*exp(2.*(double(I-1)*grid - 1.));  % estimate of wind speed in canopy using Cionco Model 
         % turb.ws(I,1)=turb.mw(I,1) *2;

         % should take the derivative of the Cionco model if we are to get
         % shear as du/dz

        
         % u/u* = uh/u* exp(a (z/H -1))

         %du/dz =a/H exp(a(z/H -1)

         turb.du(I,1)=(2/nlyr)*UH*exp(2.*(double(I-1)*grid - 1.));


         end
 
    
       figure(111)
      clf
      plot(turb.u,1:can.sze3)
 
 %        Initialize profiles for <w'u'>, <u'u'>, <v'v'>, <w'w'> and <qq>

% bottom boundary condition

      turb.uw(1,1)=0.0; 
      turb.u2(1,1)=0.;
      turb.v2(1,1)=0.;
      turb.w2(1,1)=0.;
      turb.qq(1,1)=0.;
      turb.wqq(1,1)=0.;
      turb.tau(1,1)=0.;
	  fact.fv=0.;
      
		  
	  % initialize turbulence statistics with the canopy
	  % scaled with the profile of <w'u'>

      %  seems code is not updating the uw and w2 profiles from iteration
      %  to iteration

		for I=2:(nlyr-1)
		
        turb.uw(I,1)=turb.uw(I-1,1)-1./double(nm1/3); % scales linearly from 0 to -1, ground to h
        turb.u2(I,1)=-4.*turb.uw(I,1);   % sigma_u^2 ~ 4 u* u*
        turb.v2(I,1)=-2.5*turb.uw(I,1);  % sigma_v^2 ~ 2.5 u* u*
        turb.w2(I,1)=-1.6*turb.uw(I,1);  % sigma_w^2 ~ 1.6 u* u*
        turb.qq(I,1)=turb.u2(I,1)+turb.v2(I,1)+turb.w2(I,1); % tke
        turb.q(I,1)=sqrt(turb.qq(I,1));   % sqrt(tke)
        end


		% turbulence statistics in the constant flux layer above the canopy.
		% Values scaled by <xx>/(u* u*) = constant

        for I=nlyr:can.sze3
		
        turb.uw(I,1)=-1.;  % <w'u'>/u*^2
        turb.w2(I,1)=1.6;  % <w'w'>/u*^2
        turb.v2(I,1)=2.50; % <v'v'>/u*^2
        turb.u2(I,1)=4.0;  % <u'u'>/u*^2 
        turb.qq(I,1)=8.1;  % <q'q'>/u*^2
        turb.q(I,1)=sqrt(turb.qq(I,1)); 
		end
 
		% bottom boundary conditions

      turb.qq(1,1)=0.0; 
      turb.w3(1,1)=0.;
      fact.x(1,1)=0.;
 
 
	 % initialize profiles for third order moments

      for I =2:nm1
	  turb.dw2(I,1)=(turb.w2(I+1,1)-turb.w2(I-1,1))/gridp2;
      end



  
%   Start of the main computational loop


      mxiter= 500;      % maximum number of iterations
      temp1=0; 
      temp2=0.;
      
	  turb.wwu(1,1)=0. ;
      turb.w31(1,1)=0.;
      turb.wvv1(1,1)=0.;
      turb.wuu1(1,1) = 0;

      iters=0; 
	  J=1;

     t_test=10; % initial
      pval=1;

     figure(100)
     hold on

% xlabel('xx m2 s-2')
% ylim([0 300])
% title ('turbulence variance profiles')

xlabel('u m s-1')
ylim([0 300])
title ('turbulence velocity profiles')


h=1;
 
 % Iteration loop

      
       % this works better. need to get past first few runs as the change
       % is small and run ends before it is ready.  found wt =0.5 was too
       % large and ended up leading to oscillations.  I am now just solving
       % for changes in the canopy as above the canopy there is little
       % change in iterations as well mixed.
       
     while (abs(t_test) > 2.33 && iters < mxiter || iters < 15)   %Student t value of 0.01, 2.33 one sided alpha, difference within 99%


%     Weighting coefficients for the integral solutions; 
%	  Tilden changed them during the iteration 
%     With a faster computer I am making them more strict

% I also need to question these weighting factors
% If I am using them for solving a non-linear differential equation, I want to use
% an implicit solution, that is stiffer. In this case the weighting factors approach 1.



        % weights for solving DiffEq, from Meyers, with modification

	  
         if( J <= 10)
		 wt=.01 * J;
         onemwt=1.-wt;
         else
         wt=0.1 + (J-10)*0.01;
            if wt> 0.25
             wt=0.25;
            end
          onemwt=1.-wt;
          end
        
 

        % did a bunch of plot tests. With wt near zero, I only get profiles
        % like initial conditions. With wt near one things tend to blow up
        % at the top of the domain.  wt = 0.25 seems a sweet spot and there
        % is ample weighting of the LAI profile on dwu that gets
        % incorporated into the profiles of uw, w2. u2 and v2

      
		
 
% The zero plane displacement is updated every 2 interations
 
% And so is the constant for the dissipation length scale, gamma 
 
% need to assess modulo for every other point

      if(mod(J,2) == 0)
     [can] = ZPD(turb,can); 
     [can]= LGEN(can);   % l and ld are computed here
      can.gamma=23.053*.4*(3.-can.d)/can.l(can.sze3,1);
  
	  end

	  
%  compute products that change during each interation
%  but are used several times in the loop 

 % turbulence time scale
      
    for I=2:can.sze3
	  
         turb.tau(I,1)=can.ld(I,1)/turb.q(I,1);
         fact.x(I,1)=turb.tau(I,1)/(2.*C8); 
    end
  

%  the flux divergences of w'u' is balanced by the drag force term, Cd a u^2

% d<w'u'>/dz=-Cd a(z) <u>^2  

% in this case u is normalized by u* so we have
	
%	d<w'u'/u*^2>/dz=-Cd a(z) {<u>/u*}^2

% but somewhere we lose the normalization and I am not sure where yet.

      for I=2:nm1
	  
         fact.dx(I,1)=(fact.x(I+1,1)-fact.x(I-1,1))/gridp2 ;
         turb.duw(I,1)=-can.pad(I,1)*turb.u(I,1)*turb.u(I,1)*can.drag ;
	  end

      turb.duw(1,1)=0.0;
 

 
%   Integration of the momentum equation to obtain
%   a profile of Reynolds stress <u'w'>


% At the boundary a quartic is fit through the derivatives (dfdz, f')
% and integration is performed through the first and second points


	  %f(n+1) =f(n) + dz/12 (5f'(n-1) + 8f'(n) - f'(n+1))
 
      % dy/dt, Eulers method y1 = y0 + yo(t1-t0) 
       % turb.uw(2,1)=turb.uw(2,1)+wt*(turb.uw(1,1)+grid/12.*(5.*turb.duw(1,1)+ 8.* turb.duw(2,1) ...
       %          -turb.duw(3,1))); 

      turb.uw(2,1)=onemwt*turb.uw(2,1)+wt*(turb.uw(1,1)+grid/12.*(5.*turb.duw(1,1)+ 8.* turb.duw(2,1) ...
                -turb.duw(3,1))); 

     
 
   
% The integration of the differential equation is performed using a method
% that fits a cubic polynomial through 4 derivative points (df/dz, or f').
% Integration is performed on points 2 and 3.

% The numerical formula is:

	  %f(n+1) =f(n) + dz/24 (-f'(n-1) + 13f'(n) + 13f'(n-1) - f'(n+2))

% There is also a recursive filter applied to the numerical scheme as scene
% with the weighting of f(n) and the rest of the algorithm

% Here uw is normalized by u*^2 and equals -1 for the first iteration.

% Is momentum transfer adjusting to the changes in canopy structure, as the differnce
% in canopy roughness affects shear and u*=du/dz k/(z-d)?

 
	   %f(n+1) =f(n) + dz/24 (-f'(n-1) + 13f'(n) + 13f'(n-1) - f'(n+2))

	  for I=3:(nlyr + 1)
         
            turb.uw(I,1)=onemwt*turb.uw(I,1)+wt*(turb.uw(I-1,1)+grid/24.*(-turb.duw(I-2,1)+ ...
                13.*turb.duw(I-1,1)+13.*turb.duw(I,1)-turb.duw(I+1,1)));
      end
 
 
% a constant flux layer occurs for w'u' above the canopy because d<w'u'>/dz=0

% The Reynolds shear stress at the canopy-atmosphere interface is equal to the integral
% of d<w'u'> from 0 to h.

       for(I=nlyr+2:can.sze3)
       turb.uw(I,1)=turb.uw(nlyr+1,1); 
       end

% Because the closure assumption for the pressure term produces third orderm 
% moments in the form of <wuu> C8/tau or <wwu> C8/tau by assuming steady state
% conditions d<>/dt=0 we can solve for these equations directly
	   

% <w'u'u'> = -tau/C8 (-2 <w'w'u'> d<u>/dz - <w'w'> d<u'u'>/dz -2<w'u'>d<w'u'>/dz
%	   + 2 <w'u'> <U><U> a Cd  
 
%  The third order term, wuu, not containing 
%  a gradient of uu is grouped into wuu1

      
	  for I=2:nm1
	  
         % maybe error in old form. Put in new version of equation from
         % Meyers 2025 Fortran

         % turb.wuu1(I,1)=-fact.x(I,1)*(2.*turb.uw(I,1)*turb.duw(I,1)+turb.du(I,1)*turb.wwu(I,1) ...
         % -2.*turb.uw(I,1)*turb.u(I,1)*turb.u(I,1)*can.pad(I,1)*can.drag);

            turb.wuu1(I,1)=-fact.x(I,1)*(2.*turb.uw(I,1)*turb.duw(I,1)+turb.du(I,1)*turb.wwu(I,1)...
         *(2.-0.67 *ALPH1-2*ALPH2)-2*turb.uw(I,1)*turb.u(I,1)*turb.u(I,1)*can.pad(I,1)*can.drag);
	  end
 
      turb.wuu1(can.sze3,1)=(2.*turb.wuu1(can.sze3-3,1)-9.*turb.wuu1(can.sze3-2,1)+18.*turb.wuu1(can.sze3-1,1))/11.; 


%       wwu is the transport of wu 
%       dwwu is dwwu/dz
 


% <w'w'u'> = -tau/C8 (- <w'w'w'> d<u>/dz - <w'u'> d<w'w'>/dz - <w'w'> Cd a <u><u> 
%  -2<w'w'>d<w'u'>/dz


       for I=2:nm1
	   
         turb.wwu(I,1)=turb.wwu(I,1)*onemwt-wt*fact.x(I,1)*(2.*turb.w2(I,1)*turb.duw(I,1) ...
               +turb.uw(I,1)*turb.dw2(I,1)-turb.w2(I,1)*turb.u(I,1)*turb.u(I,1)*can.pad(I,1)*can.drag  ...
                +turb.du(I,1)*turb.w3(I,1)*(1-ALPH1-ALPH2)-(ALPH0-0.67*ALPH1)...
                *turb.wqq(I,1)*turb.du(I,1)-(ALPH1-ALPH2)*turb.wuu(I,1)*turb.du(I,1));

        turb.w31(I,1)=-fact.x(I,1)*((-0.67*ALPH1+2*ALPH2)*turb.wwu(I,1)*turb.du(I,1));
        turb.wvv1(I,1)=-fact.x(I,1)*(1.33*ALPH1*turb.wwu(I,1)*turb.du(I,1));

	   end

 
        turb.wwu(can.sze3,1)=(2.*turb.wwu(can.sze3-3,1)-9.*turb.wwu(can.sze3-2,1)+18.*turb.wwu(can.sze3-1,1))/11.;
        turb.w31(can.sze3,1)=(2.*turb.w31(can.sze3-3,1)-9.*turb.w31(can.sze3-2,1)+18.*turb.w31(can.sze3-1,1))/11.;
        turb.wvv1(can.sze3,1)=(2.*turb.wvv1(can.sze3-3,1)-9.*turb.wvv1(can.sze3-2,1)+18.*turb.wvv1(can.sze3-1,1))/11.;


% d <w'w'u'>/dz


       for I=2:nm1
        turb.dwwu(I,1)=(turb.wwu(I+1,1)-turb.wwu(I-1,1))/gridp2;
         turb.dw31(I,1)=(turb.w31(I+1,1)-turb.w31(I-1,1))/gridp2;
          turb.dwvv1(I,1)=(turb.wvv1(I+1,1)-turb.wvv1(I-1,1))/gridp2;
       end
 
       turb.dwwu(can.sze3,1)=0.; 
       turb.dwwu(2,1)=0.;   % in Meyers new code



% solving for du/dz from algebraic manipulation of d<u'w'>/dt budget

       for I=2:can.sze3
	   
          fact.denom(I,1)=ALPH0*turb.qq(I,1)-turb.w2(I,1)+ALPH1*(turb.w2(I,1)-turb.qq(I,1)/3.) ...
         +(ALPH1-ALPH2)*(turb.u2(I,1)-turb.qq(I,1)/3.)+ALPH2*(turb.w2(I,1)-turb.qq(I,1)/3.);

         % do not divide by zero
 
		  if (abs(fact.denom(I,1)) < 1e-6)
		  fact.denom(I,1)=1.; 
          %fact.denom(I,1)=1e-6; 
          end
       
		  turb.du(I,1)=(turb.dwwu(I,1)+C3*turb.uw(I,1)/turb.tau(I,1))/fact.denom(I,1);
        end
	   
          
    
  %   Constants used several times in the loop 

	   for I=2:nm1
	   
		 can.ld(I,1)=can.l(I,1)*can.gamma;
         fact.x1(I,1)=fact.x(I,1)*turb.w2(I,1)/gridsq;
         fact.x2(I,1)=(fact.x(I,1)*turb.dw2(I,1)+fact.dx(I,1)*turb.w2(I,1))/gridp2 ;
         fact.x3(I,1)=C2/turb.tau(I,1);
         fact.x4(I,1)=turb.qq(I,1)*C2/(3.*turb.tau(I,1)) ;
         fact.x5(I,1)=2.*power(turb.q(I,1),3.)/(3.*can.ld(I,1)) ;
         fact.x6(I,1)=turb.uw(I,1)*turb.du(I,1);
         fact.x7(I,1)=2.*power(turb.u(I,1),3.)*can.pad(I,1)*can.drag;
         turb.dwuu1(I,1)=(turb.wuu1(I+1,1)-turb.wuu1(I-1,1))/gridp2;
         turb.dwqq(I,1)=-15./9.*(2.*fact.x6(I,1)-fact.x7(I,1)+3.*fact.x5(I,1));
		 end


      can.ld(can.sze3,1)=can.l(can.sze3,1)*can.gamma; 


	  % Variance terms are solved using central difference scheme for derivatives
	  % and using boundary problem solution methods

	  % Non-linearities arise when substituting third order terms into second order budget eq


% Compute the U2 budget Profile


	  % d<u'u'>/dt = (-2 <u'u'> d<u>/dz - d<w'u'u'>/dz + 2 <U><U><U> a Cd  + ?????

		      
	  for I=2:nm1
	  
         thom.a1(I,1)=fact.x1(I,1)-fact.x2(I,1);
         thom.b1(I,1)=-2.*fact.x1(I,1)-fact.x3(I,1);
         thom.c1(I,1)=fact.x1(I,1)+fact.x2(I,1);
         thom.d1(I,1)=fact.x6(I,1)*(2.-2.*ALPH2-.67*ALPH1) -fact.x4(I,1)+fact.x5(I,1)- ...
            2./15.*turb.dwqq(I,1) -fact.x7(I,1) +turb.dwuu1(I,1); 

         % found typo in old matlab, new meyers code has dwuu1.. i was not using that
         % term which raised my curiosity. effect is tiny
         % old line
           
	  end
 
    
% Call Thomas algorithm for simultaneous Eq
	  
     [thom] =THOMAS(thom,can);

      ratio=-4.0*turb.uw(can.sze3,1)/turb.u2(can.sze3,1);
      turb.u2(can.sze3,1)=-4.0*turb.uw(can.sze3,1); % here the absolute value of <u'u'> is computed

	  

    for I=1:nm1
	     K=can.sze3-I; 
         turb.u2(K,1)=abs(onemwt*turb.u2(K,1)*ratio+wt*(thom.delta(K,1)-thom.f(K,1)*turb.u2(K+1,1))); 
	end
 


%     Compute V2 Budget 

% +turb.dwvv1(I,1); was missing. in new code

	
	for I=2:nm1
         thom.d1(I,1)=-fact.x4(I,1)+1.33*ALPH1*fact.x6(I,1)-2./15 *turb.dwqq(I,1)+fact.x5(I,1)...
             +turb.dwvv1(I,1); 
    end
 
      [thom]=THOMAS(thom,can);

      ratio=-2.50*turb.uw(can.sze3,1)/turb.v2(can.sze3,1); 
      turb.v2(can.sze3,1)=-2.50*turb.uw(can.sze3,1);  % here the absolute value of <v'v'> is computed
		

      for I=1:nm1
	  
         K=can.sze3-I; 
         turb.v2(K,1)=abs(onemwt*turb.v2(K,1)*ratio+wt*(thom.delta(K,1)-thom.f(K,1)*turb.v2(K+1,1))); 
	  end
 
 

%     Compute the ww and qq budgets 

% +turb.w31(I,1) was missing

      
	  for I=2:nm1
	  
         thom.a1(I,1)=3.*thom.a1(I,1); 
         thom.b1(I,1)=-6.*fact.x1(I,1)-fact.x3(I,1);
         thom.c1(I,1)=3.*thom.c1(I,1) ;
         thom.d1(I,1)=-fact.x4(I,1)-2./15. *turb.dwqq(I,1)+fact.x5(I,1)...
              +fact.x6(I,1)*(-.67*ALPH1+2.*ALPH2) +turb.w31(I,1); 
	  end
 

      [thom]=THOMAS(thom, can);

      ratio=-1.6*turb.uw(can.sze3,1)/turb.w2(can.sze3,1);
      turb.w2(can.sze3,1)=-1.6*turb.uw(can.sze3,1); % here the absolute value of <w'w'> is computed
     	  
	  turb.qq(can.sze3,1)=turb.u2(can.sze3,1)+turb.v2(can.sze3,1)+turb.w2(can.sze3,1); 
      turb.q(can.sze3,1)=sqrt(turb.qq(can.sze3,1)); 
 
      for I=1:nm1
	  
         K=can.sze3-I; 
         turb.w2(K,1)=abs(onemwt*turb.w2(K,1)*ratio+wt*(thom.delta(K,1)-thom.f(K,1)*turb.w2(K+1,1))); 
         turb.qq(K,1)=turb.u2(K,1)+turb.v2(K,1)+turb.w2(K,1);
         turb.q(K,1)=sqrt(turb.qq(K,1));
	  end



% Looking at first principles these terms are functions of them selves

	  % better to solve with implicit differential equations and re-pose these terms



        for I=2:nm1
     	   

 % d<w'w'>/dz
          turb.dw2(I,1)=(turb.w2(I+1,1)-turb.w2(I-1,1))/gridp2;

 % d<u'u'>/dz
          turb.du2(I,1)=(turb.u2(I+1,1)-turb.u2(I-1,1))/gridp2; 

 % d<v'v'>/dz         
		  turb.dv2(I,1)=(turb.v2(I+1,1)-turb.v2(I-1,1))/gridp2; 

 % <w'u'u'>=-tau/C8(2 <w'w'u'> d<u>/dz+ 2<w'u'> d<w'u'>/dz +2 <u'w'> Cd a(z) <u>^2  

		  turb.wuu(I,1)=-fact.x(I,1)*(turb.w2(I,1)*turb.du2(I,1)+2.*turb.uw(I,1)*turb.duw(I,1)...
          +turb.du(I,1)*turb.wwu(I,1)...
         *(2-0.67*ALPH1-2*ALPH2)...
         -2.*turb.uw(I,1)*turb.u(I,1)*turb.u(I,1)*can.pad(I,1)*can.drag)*wt ...
         +onemwt*turb.wwu(I,1);
          
% <w'v'v'>=-tau/C8(<w'w'> d<v'v'>/dz

		  % d <w'v'v'>/ dt ~ -tau/C8(<w'w'> d<v'v'>/dz


		 turb.wvv(I,1)=-fact.x(I,1)*(turb.w2(I,1)*turb.dv2(I,1) +1.33*ALPH1*turb.wwu(I,1)*turb.du(I,1))*wt...
             +onemwt*turb.wvv(I,1);

            
         
% <w'w'w'>=-tau/C8(3<w'w'> d<w'w'>/dz

		 turb.w3(I,1)=-fact.x(I,1)*(3.*turb.w2(I,1)*turb.dw2(I,1) +(-0.67*ALPH1 +2*ALPH2)*turb.wwu(I,1)*turb.du(I,1))*wt...
             + onemwt*turb.w3(I,1); 

    
% <w'w'u'> + <w'v'v'> + <w'w'w'>
         
		 turb.wqq(I,1)=turb.wuu(I,1)+turb.wvv(I,1)+turb.w3(I,1);
       end
 

% **************************************************************
%      Determine the third order moments at the top of the grid
%      assuming the gradient of the triple moment is 0
% **************************************************************



	    turb.w3(can.sze3,1)=(2.*turb.w3(can.sze3-3,1)-9.*turb.w3(can.sze3-2,1)+18.*turb.w3(can.sze3-1,1))/11.;
        
	    turb.wvv(can.sze3,1)=(2.*turb.wvv(can.sze3-3,1)-9.*turb.wvv(can.sze3-2,1)+18.*turb.wvv(can.sze3-1,1))/11.;
        
	    turb.wuu(can.sze3,1)=(2.*turb.wuu(can.sze3-3,1)-9.*turb.wuu(can.sze3-2,1)+18.*turb.wuu(can.sze3-1,1))/11.;
        
	    turb.wqq(can.sze3,1)=turb.w3(can.sze3,1)+turb.wvv(can.sze3,1)+turb.wuu(can.sze3,1); 


 
% The mean horizontal wind profile (mw) is computed by
% integrating from the top of the domain downward du, which should not be normalized by u*
% as it is computed from the budget equations
	   
% a logarithmic wind profile is assumed in the lowest layers of the canopy
	   

 turb.u(nm1,1)=onemwt*turb.u(nm1,1)+wt*(turb.u(can.sze3,1)-grid/12.*(5.*turb.du(can.sze3,1)+8.*turb.du(nm1,1)- ...
     turb.du(nm2,1)));
 

      if(turb.u(nm1,1) < 0.)
		  turb.u(nm1,1)=0.;
      end


	  % The numerical integration formula is:

	  %f(n+1) =f(n) + dz/24 (-f'(n-1) + 13f'(n) + 13f'(n-1) - f'(n+2))

 
      % from top to bottom in this loop

        



          for I=2:nm3

          K=can.sze3-I ;

           % Adams Bashford Integration, requires information from 3 levels 
          %f(n+1) =f(n) + dz/24 (-f'(n-1) + 13f'(n) + 13f'(n-1) - f'(n+2))

          % turb.u(K,1)=onemwt*turb.u(K,1)+wt*(turb.u(K+1,1)-grid/24.*(-turb.du(K+2,1)+13.*turb.du(K+1,1) ...
          % +13.*turb.du(K,1)-turb.du(K-1,1)));

% revised with numerics

           turb.u(K+1,1)=onemwt*turb.u(K+1,1)+wt*(turb.u(K,1)+grid/24.*(-turb.du(K-1,1)+13.*turb.du(K,1) ...
           +13.*turb.du(K-1,1)-turb.du(K+2,1)));


         if(turb.u(K,1) <0.)
			 turb.u(K,1)=0.;
         end

         end

  	 
         % try simpler numeric for lower boundary at n=2

         % Adams Bashford used above..but requires 3 levels
         %f(n+1) =f(n) + dz/12 (5f'(n-1) + 8f'(n) - f'(n+1))

         % try simpler integration, Trapezoid rule

         % (b-a) 0.5(f(a) + f(b))

      	
	  I=2; 
      % turb.u(I,1)=onemwt*turb.u(I,1)+wt*(turb.u(I+1,1)-grid/24.*(9.0*turb.du(I,1)+19.0 * ...
		  % 	 turb.du(I+1,1)-5.0*turb.du(I+2,1)+turb.du(I+3,1))); 

          % turb.u(I,1)=onemwt*turb.u(I,1)+wt*(turb.u(I-1,1)+grid/12.*(5.0*turb.du(I-1,1)+8.0 * ...
		   % 	 turb.du(I-1,1)-turb.du(I,1))); 

            turb.u(I,1)=onemwt*turb.u(I,1)+wt*(grid*(turb.u(I-1,1)+turb.u(I,1))/2);




     
      % 
	 if(turb.u(I,1) < 0.)
	  turb.u(I,1)=0.0 ;
     end
     
	 
      fact.fv=turb.u(2,1)*.4/(log((can.zg+grid)/can.zg)); 
      turb.uw(1,1)=onemwt*turb.uw(1,1)-wt*fact.fv*fact.fv; 


	  turb.ustar=sqrt(abs(turb.uw(nlyr+1,1)));
      
	 % try testing for convergence of the integrated tke profile
	 % using paired t test for convergences  
	  
	  [h,pval,ci,stats]= ttest(turb.old(1:can.sze1),turb.qq(1:can.sze1));

      sumdiffqq=sum((turb.old(1:can.sze1)-turb.qq(1:can.sze1)));

      t_test=stats.tstat;

	  % give new to old
      
	  turb.old=turb.qq;
      
	 	  %  check for convergence  
      
%plot(turb.uw(1:300),1:300,'.-','MarkerSize', 10)
%plot(turb.qq(1:300),1:300,'.-','MarkerSize', 10)

plot(turb.u(1:300),1:300,'.-','MarkerSize', 10)



		
		iters =iters +1; 
		J = J +1;
	 
 	end
        % } while (fabs(ttest) > 2.33);  %P 0.01

	 
 
% call zero plane displacement

     [can]= ZPD(turb,can);
 
      

 
     
% z0/h
% note turb.u is u/u*

	  
 % normalized by momentum stress, u*

 % Note, if I normalize the higher order moments by u* then I match
 % the upper boundary conditions, except for u

% logic may be circular if u* is 1.0 anyhow

      for I=1:can.sze3
	  				  
         turb.us(I,1)=turb.u(I,1)*(1.+turb.v2(I,1)/(2.*turb.u(I,1)*turb.u(I,1)));
		 turb.v(I,1)=power((turb.us(I,1)*turb.us(I,1)-turb.u(I,1)*turb.u(I,1)),0.5);
		 turb.u(I,1)=turb.u(I,1)/turb.ustar; 
         turb.u2(I,1)=turb.u2(I,1)/power(turb.ustar,2);
         turb.v2(I,1)=turb.v2(I,1)/power(turb.ustar,2);
         turb.w2(I,1)=turb.w2(I,1)/power(turb.ustar,2);
         turb.uw(I,1)=turb.uw(I,1)/power(turb.ustar,2);
         turb.du(I,1)=turb.du(I,1)/turb.ustar;
         turb.q(I,1)=turb.q(I,1)/turb.ustar;
		 turb.duw(I,1)=turb.duw(I,1)/power(turb.ustar,2);
         turb.w3(I,1)=turb.w3(I,1)/power(turb.ustar,3);
         turb.wuu(I,1)=turb.wuu(I,1)/power(turb.ustar,3);
         turb.wvv(I,1)=turb.wvv(I,1)/power(turb.ustar,3);
         turb.wwu(I,1)=turb.wwu(I,1)/power(turb.ustar,3);

         turb.qq(I,1)=turb.q(I,1)*turb.q(I,1);
         turb.tau(I,1)=turb.qq(I,1)*can.ld(I,1)/power(turb.qq(I,1),1.5);
         fact.x(I,1)=turb.tau(I,1)/(2.*C8);


	  end



%????? why is can.z0 computed as a profile. check c code

% doesnt seem to be used so comment out on debug 2/5/24


  % for I=can.sze1:can.sze3	
  % %  can.z0(I,1)=(double(I-1)*grid-can.d)/exp(0.4*turb.u(I,1));
  % can.z0(I,1)=(double(I-1)*grid-can.d)/exp(0.4*turb.u(I,1)/turb.ustar);
  % % Z0(I)=(FLOAT(I-1)*G-DTHOM)/EXP(MW(I)*.4/RUST)   
  %  end

	  % dissipation =qq^(3/2)/aL

	  % turbulence time scale  = qq/epsilon
 
   
       turb.us(1,1)=0.;
       turb.u(1,1)=0.;
       turb.du(1,1)=0.0;
      
      % plot the profiles

      plotprofiles(turb, can);

      toc;

	  % end of main

