
function [asave,bsave] = RungeKuttaFuture_schwartz3d( nrunge ,Trunge, A, B, omega0  )
        
    h = Trunge/nrunge;
    
    dya = @(t,b)  b*A + 0.5 * b*omega0*b';
    dyb = @(t,b) [b*B(:,1), b*B(:,2), b*B(:,3) ];
    
    %Initial conditions
    
    a = 0; 
    b = [1 0 0];
    asave = a;
    bsave = b;
    
    %Fourth order
    for t = 0+h : h : Trunge %parte da 0+h perchè ho già valore in zero
        
        a1 = dya( t , b );
        b1 = dyb( t , b );
    
        a2 = dya( t+h/2, b+h/2*a1 );
        b2 = dyb( t+h/2, b+h/2*b1 );
    
        a3 = dya(t+h/2, b+h/2*a2 );
        b3 = dyb(t+h/2, b+h/2*b2 );
        
        a4 = dya(t+h,   b+h*a3 );
        b4 = dyb(t+h,   b+h*b3 );
    
        a = a + h/6 * (a1 + 2*a2 + 2*a3 + a4);
        b = b + h/6 * (b1 + 2*b2 + 2*b3 + b4);
        asave = [asave; a];
        bsave = [bsave; b];
    
    end 
    asave = asave(end);
    bsave = bsave(end,:);
end



%Exact solutions of a and b are in Eq. 27-28 Schwartz

% A = (kd*mud + sigx*sigd*rhoxd)*( (1-exp(-kd*T)) - kd*T  )/(kd^2);
% 
% B = sigd^2 * ( 4*(1-exp(-kd*T)) - (1-exp(-2*kd*T)) -2*kd*T )/(4*kd^3);
% 
% C = ((kr*mur + sigx*sigr*rhoxr)*( (1-exp(-kr*T)) - kr*T))/kr^2;
% 
% D = sigr^2 * ( 4*(1-exp(-kr*T)) - (1-exp(-2*kr*T)) - 2*kr*T)/(4*kr^3);
% 
% E = (  (1-exp(-kd*T))+ (1-exp(-kr*T)) - (1-exp(-(kd+kr)*T)) )/ (kr*kd*(kr+kd));
% 
% F = ( kd^2*(1-exp(-kr*T)) + kr^2 * (1-exp(-kd*T)) - kd*kr^2*T - kr*kd^2*T) / (kr^2*kd^2*(kr+kd));
% 

%a = A-B-C-D+sigd*sigr*rhoxr*(E+F)
%Z = [-(1-exp(-kd*T))/kd,  (1-exp(-kr*T))/kr]
