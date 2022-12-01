function D = pearson_cc(X,Y)
        
    X=X';
    Y=Y';
    N=length(X);
    XY=X.*Y;
    X2 = X.*X;
    Y2=Y.*Y;
    num = (N.*(sum(XY)))-(sum(X).*sum(Y));
    densq = (N.*sum(X2)-(sum(X).^2)).*(N.*sum(Y2)-(sum(Y).^2));
    den = sqrt(densq);
    D = num./den;
    
end