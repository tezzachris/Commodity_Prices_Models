
%constant risk premia
function [ll,lp, Xt, et, a, Z] = loglik_schwartz2d( theta, yt, TrungeF, nrunge, h)
    
    n = size(yt,1);  %n.obs
    ny = size(yt,2); %n.futures

    %load the paramters
    theta_x = num2cell(theta);
    [mux,r,... 
     mud,kd, ...
     sigx,rhoxd,sigd,... 
     lambda] = deal(theta_x{1:end-ny});
    sigepsi = theta(end-ny+1 : end);
    H = diag( sigepsi.^2 );
   
    %risk neutral
    A = [ r - sigx^2/2 ; 
          mud * kd ];

    %physical
    Af = [mux - sigx^2/2 ; 
          mud * kd - lambda ];
     
    B = [ 0 , -1;
          0 , -kd];

    Bf = B;
            
    omega0 =  [ sigx^2     rhoxd*sigx*sigd  ;
                rhoxd*sigx*sigd    sigd^2   ];
    
    %Runge kutta per coefficienti affini
    a = []; Z = [];
    for i = 1:length(TrungeF)
         [ait,bit] = RungeKuttaFuture_schwartz2d(nrunge, TrungeF(i), A, B, omega0);
         a = [a, ait]; 
         Z = [Z; bit];
    end
    
    %Exact solution is  
    % alphah = mud - lambda/kd;
    % 
    % a = ( r - alphah + 1/2 * sigd^2/kd^2 - (sigx*sigd*rhoxd)/kd )  * TrungeF ...
    %                   + 1/4*sigd^2 * ((1-exp(-2*kd*TrungeF))/(kd^3)) +...
    %                    ( alphah*kd + (sigx*sigd*rhoxd) -  sigd^2/kd ) * ((1-exp(-kd*TrungeF))/(kd^2)) ;
    % 
    % Z = [1 , -(1-exp(-kd*TrungeF))/kd] ;

    d = size(A,1); %number of state variables
    ll = 0; 
    lp = zeros(n,1);
    Xt = zeros(n,d); 
    et = zeros(n,ny); %error terms 
    I = eye(d);
    
    %Initial value for state variables
    Xt0(1) = Af(1) ;
    Xt0(2) = Af(2)/kd ;
    Pt0(1,1) = sigx^2 / 2; 
    Pt0(2,2) = sigd^2 / (2*kd); 
    Pt0(1,2) = rhoxd*sigd*sigx ; 
    Pt0(2,1) = Pt0(1,2);
    Pt = Pt0;

    for i = 1 : n 
        
            if i > 1
                Xt0 = Xt(i-1,:);
            end

            %Predicted state
            Xpred = h * Af' + ( (I + h * Bf) * Xt0' )'  ; %mean
            Ppred = (I+h*Bf)* Pt *(I+h*Bf)' + h*omega0; %variance

            %Predicted obs
            ypred = a + ( Z * Xpred' )' ; %the risk neutral drift is observable
            
            %Prediction error
            vt = yt(i,:) - ypred ; %mean
            Ft = Z * Ppred * Z' + H; %variance
            

            % [~,d] = eig(Ft); d = min(d);
            % 
            % if d(1) ~= 0
            %      'problema invertibilità'
            % end

            Kt = (Ppred * Z')/(Ft);
            Xt(i,:) = Xpred + (Kt * vt')'; 
            Pt =  Ppred - Kt*Z*Ppred;
            
            lp(i) = 0.5 * ( ny * log(2*pi) + log( det(Ft) ) +  (vt / Ft * vt') ); %vt needs to be column vector
            ll = ll + lp(i);  
            et(i,:) = vt/sqrt(Ft);
    end   
end