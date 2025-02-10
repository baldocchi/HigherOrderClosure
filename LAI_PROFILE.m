function [can] = LAI_PROFILE(can)
   

           
          can.ht=10.;			% canopy height, m
		  can.lai=4.0;	        % leaf area index

% need to modify for oak savanna and put out normalzed LAI profile
				  

            % Evaluate how LAI and other canopy structural variables vary
            % with time

              
         jtot=can.sze1-1;

        % height of mid point of layer scaled of normalized height

        ht_midpt(1) = .05;
        ht_midpt(2)= .15;
        ht_midpt(3)= .25;
        ht_midpt(4) = .35;
        ht_midpt(5) = .45;
		ht_midpt(6) = .55;
        ht_midpt(7)= .65;
        ht_midpt(8)= .75;
        ht_midpt(9) = .85;
        ht_midpt(10) = .95;

 

      % lai of the layers at the midpoint of height


        lai_freq(1) = .025;
        lai_freq(2) = .025;
        lai_freq(3) = .05;
        lai_freq(4) = .1;
        lai_freq(5) = .2;
        lai_freq(6) = .3;
        lai_freq(7) = .2;
        lai_freq(8) = .05;
        lai_freq(9) = .025;
        lai_freq(10) = .025;
        

		% lai_freq(1) = .025;
        % lai_freq(2) = .025;
        % lai_freq(3) = .025;
        % lai_freq(4) = .025;
		% lai_freq(5) = .05;
        % lai_freq(6) = .1;
        % lai_freq(7) = .2;
        % lai_freq(8) = .3;
        % lai_freq(9) = .2;
        % lai_freq(10) = .05;
        % 
     

% /*
%    Beta distribution
% 
%    f(x) = x^(p-1) (1-x)^(q-1) / B(v,w)
% 
%    B(v,w) = int from 0 to 1 x^(p-1) (1-x)^(q-1) dx
% 
%   p = mean{[mean(1-mean)/var]-1}
% 
%   q =(1-mean){[mean(1-mean)/var]-1}
% 
%   *****************************************************************
% */
        TF = 0.;
        MU1 = 0.;
        MU2 = 0.;
        integr_beta = 0.;


%  Height at the midpoint


        for I = 1:10
            
        % Normalize height 
     
       % ht_midpt(I) /= can.ht;   

        
        % Total F in each layer. Should sum to LAI
        

        TF = TF + lai_freq(I);

        
        % weighted mean lai
        
        MU1 = MU1 + (ht_midpt(I) * lai_freq(I));

        
        % weighted variance
        
        MU2 = MU2 +  (ht_midpt(I)* ht_midpt(I) * lai_freq(I));
        end  % next I 

        
        % normalize mu by lai
        
        MU1 =MU1/ TF;
        MU2 =  MU2/TF;

        
        % compute Beta parameters
        

        P_beta = MU1 * (MU1 - MU2) / (MU2 - MU1 * MU1);
        Q_beta = (1. - MU1) * (MU1 - MU2) / (MU2 - MU1 * MU1);
        P_beta = P_beta -1;
        Q_beta = Q_beta - 1.;

        % /*
        % '  integrate Beta function, with Simpson's Approx.
        % '
        % '  The boundary conditions are level 1 is height of ground
        % '  and level jtot+1 is height of canopy.  Layer 1 is between
        % '  height levels 1 and 2.  Layer jtot is between levels
        % '  jtot and jtot+1
        % '
        % '  Thickness of layer
        % */

        dx = 1. / jtot;

        DX2 = dx / 2.;
        DX4 = dx / 4.;
        X = DX4;

        F2 = (power(X,P_beta)) *power((1. - X), Q_beta);
        X = X + DX4;
        F3 = power(X, P_beta) *power((1. - X),Q_beta);

        
        % start integration at lowest boundary
        

        beta_fnc(1) = DX4 * (4. * F2 + F3) / 3.;
        integr_beta = integr_beta+ beta_fnc(1);

        JM1=jtot-1;

        for I = 2:JM1
        
        F1 = F3;
        X = X+ DX2;
        F2 = power(X, P_beta) * power((1. - X), Q_beta);
        X = X + DX2;
        F3 = power(X, P_beta) * power((1. - X), Q_beta);
        beta_fnc(I) = DX2 * (F1 + 4. * F2 + F3) / 3.;
        integr_beta =integr_beta + beta_fnc(I);
        end

        F1 = F3;
        X = X + DX4;

        F2 = power(X, P_beta) * power((1. - X),Q_beta);

        
        %  compute integrand at highest boundary
        

        beta_fnc(jtot) = DX4 * (F1 + 4. * F2) / 3.;

        integr_beta = integr_beta + beta_fnc(jtot);


% /*
%         '   lai_z IS THE LEAF AREA AS A FUNCTION OF Z
%         '
%         '   beta_fnc is the pdf for the interval dx
% */

        can.lai_z(1,1) = beta_fnc(1) *can.lai / integr_beta;

        for I = 2:JM1
        can.lai_z(I,1) = beta_fnc(I) *can.lai / integr_beta;
        end

        can.lai_z(jtot,1) = beta_fnc(jtot) * can.lai / integr_beta;

        cum_ht = 0;
        cum_lai = 0;

  
        for I = 1:jtot
        
        % /*
        % ' re-index layers of lai_z.
        % ' layer 1 is between ground and 1st level
        % ' layer jtot is between level jtot and top of canopy (jtot+1)
        % */
        % 
        % %cum_ht += delz;
        % %cum_lai += can.lai_z(I);
        % 
        
        dz=can.ht/jtot;   % distance per layer

        can.pad(I,1)=can.lai_z(I,1)/dz;  % lai density m2 m-3

       
        end  % next I 

           
    
end