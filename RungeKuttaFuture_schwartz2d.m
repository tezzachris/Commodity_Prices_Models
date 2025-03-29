
function [asave,bsave] = RungeKuttaFuture_schwartz2d( nrunge ,Trunge, A, B, omega0  )
       
    h = Trunge/nrunge;
    
    dya = @(t,b)  b*A + 0.5 * b*omega0*b';
    dyb = @(t,b) [b*B(:,1), b*B(:,2)];
    
    %Initial conditions
    
    a = 0; 
    b = [1 0];

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
    %salvare solo valori a convergenza
    asave = asave(end);
    bsave = bsave(end,:);
end