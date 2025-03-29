
%Log-Likelihood for Schwartz one-dimensional model for commodity log-spot prices

%Inputs

%theta: vector parameters 
%yt: matrix of log-future prices (time per row, maturity per column)
%T: maturities vector (expressed in years)
%h: Euler discretization time step

%Outputs
%ll: log-likelihood value
%lp: log-likelihood point in time values
%Xt: filtered state variable/s
%et: filtered errors
%a,Z: affine coefficients in log-price equation

function [ll, lp, Xt, et,a,Z] = loglik_schwartz1d(theta, yt, T, h)
    
    n = size(yt,1); %n observations 
    ny = size(yt,2); %n futures

    %load the paramters
    theta_x = num2cell(theta);

    [mux,r,... 
     kx,sigx] = deal(theta_x{1:end-ny});
    
    sigepsi = theta(end-ny+1 : end);
    
    H = diag( sigepsi.^2 );

    %assume constant risk premia B=Bf, that is lambda
    A = kx * r - sigx^2/2;

    Af =  kx * mux  -  sigx^2 / 2 ;
     
    B =  -kx ;
            
    Bf = B;

    omega0 = sigx^2 ;
    
    %Runge kutta for affine coefficients
    alphastar = A/kx ;
    a = (1-exp(-kx*T))*alphastar + sigx^2/(4*kx) * (1-exp(-2*kx*T));
    Z = exp(-kx*T);

    %Linear Kalman Filter equations

    ll = 0; 
    lp = zeros(n,1);
    Xt = zeros(n,1); 
    I = 1;
    et = zeros(n,ny);

    %initial state variables
    Xt0 = mux  - sigx^2 /2;
    
    Pt = sigx^2/(2*kx) * (1 - exp(-2*kx*mean(T))); 
    
    for i = 1 : n 
        if i > 1
           Xt0 = Xt(i-1,:);
        end
        %Predicted state
        Xpred = h * Af' + ( (I + h * Bf) * Xt0' )'; %E[Xt|yt-1]
        Ppred = (I+h*Bf)* Pt *(I+h*Bf)' +  h * omega0 ; %Var[Xt|yt-1]
        ypred = a + Z * Xpred ; %the risk neutral drift is observable
    
        %Prediction error
        vt = yt(i,:) - ypred ; %mean
        Ft = Z * Ppred * Z' + H; %variance + measurement error
        
        %Kalman gain
        Kt = (Ppred * Z)/(Ft);
        Xt(i,:) = Xpred + (Kt * vt'); 

        %Update eqs
        Pt =  Ppred - Kt*Z'*Ppred;

        lp(i) = 0.5 * ( ny * log(2*pi) + log( det(Ft) ) +  (vt / Ft * vt') ); %vt needs to be column vector
        ll = ll + lp(i);    

        et(i,:) = vt/sqrt(Ft); 

    end     
end   
