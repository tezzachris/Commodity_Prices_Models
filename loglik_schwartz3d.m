function [ll,lp,Xt,et,a,Z] = loglik_schwartz3d(theta, yt, TrungeF, nrunge, h)
    
    n = size(yt,1);   %n.obs
    ny = size(yt,2);  %n.futures
    
    %load the paramters
    theta_x = num2cell(theta);
    [mud,mur,... 
     kd,kr, ...
     sigx,rhoxd,rhoxr,...
     sigd,rhodr,sigr,... 
     lambda] = deal(theta_x{1:end-ny});
    
    sigepsi = theta(end-ny+1 : end);

    H = diag( sigepsi.^2 );

    %risk neutral
    A = [   - 0.5*sigx^2 ; 
              mud*kd ;
              mur*kr];
    %physical
    Af = [ -0.5*sigx^2 ; 
            mud * kd - lambda ;
            mur*kr];
     
    B = [ 0 ,  -1, + 1;
          0 , -kd,   0;
          0 ,  0 , -kr];

    Bf = B; %costant risk premia as in Schwartz
            
    omega0 = [   sigx^2  rhoxd*sigx*sigd rhoxr*sigx*sigr ;
                 rhoxd*sigx*sigd sigd^2  rhodr*sigd*sigr ;
                 rhoxr*sigx*sigr rhodr*sigd*sigr sigr^2];
    
    
    %Runge kutta per coefficienti affini per formula log-prezzi
    a = []; Z = [];

    for i = 1:length(TrungeF)
        [ai,Zi] = RungeKuttaFuture_schwartz3d(nrunge,TrungeF(i),A,B,omega0);
        a = [a, ai]; 
        Z = [Z; Zi];
    end

       
    %Kalman Filter

    d = 3; %number of state variables
    ll = 0; 
    lp = zeros(n,1);
    Xt = zeros(n,d); 
    I = eye(d,d);
    et = zeros(n,ny);

    %Initials

    Xt0 = [mur - mud - sigx^2/2, mud, mur];
    
    Pt0=zeros(3,3);

    Pt0(1,1) = sigx^2 / 2; 
    Pt0(2,2) = sigd^2 / (2*kd); 
    Pt0(3,3) = sigr^2 / (2*kr); 
    Pt0(1,2) = rhoxd*sigd*sigx ; Pt0(2,1) = Pt0(1,2);
    Pt0(1,3) = rhoxr*sigr*sigx ; Pt0(3,1) = Pt0(1,3);
    Pt0(2,3) = rhodr*sigd*sigr ; Pt0(3,2) = Pt0(2,3);
    Pt = Pt0;

    Pt = eye(d,d);
    
    for i = 1 : n 
        if i > 1
                Xt0 = Xt(i-1,:);
        end
            %Prediction equations
            Xpred = h * Af' + ( (I + h * Bf) * Xt0' )'  ; %mean
            Ppred = (I+h*Bf)* Pt *(I+h*Bf)' + h * omega0 ; %variance
            ypred = a + ( Z * Xpred' )'; %the risk neutral drift is observable
            
            %Prediction error
            vt = yt(i,:) - ypred ; %mean
            Ft = Z * Ppred * Z' + H; %variance
            
            %Update equations
            Kt = (Ppred * Z')/(Ft);
            Xt(i,:) = Xpred + (Kt * vt')'; 
            Pt =  Ppred - Kt*Z*Ppred;
            
            lp(i) = 0.5 * ( ny * log(2*pi) + log( det(Ft) ) +  (vt / Ft * vt') ); %vt needs to be column vector
            ll = ll + lp(i);
            et(i,:) = vt/sqrt(Ft); %prediction errors
                 
    end   
end