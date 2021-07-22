function [Value,pValue,ds] = getStats(X,A,B,C,D,SSE,statLabel,withinNames,betweenNames)
% Adapted from RepeatedMeasuresModel.m

% Hypothesis matrix H
% H = (A*Beta*C - D)'*inv(A*inv(X'*X)*A')*(A*Beta*C - D);
% q = rank(Z);
[H,q] = makeH(A,B,C,D,X);

% Error matrix E
E = C'*SSE*C;

p = rank(E+H);
s = min(p,q);
v = size(X,1)-rank(X);
if p^2+q^2>5
    t = sqrt( (p^2*q^2-4) / (p^2+q^2-5));
else
    t = 1;
end
u = (p*q-2)/4;
r = v - (p-q+1)/2;
m = (abs(p-q)-1)/2;
n = (v-p-1)/2;

switch statLabel
    case 'Wilks'
        % ~~~ Wilks' Lambda = L
        % Formally, L = |E| / |H+E|, but it is more convenient to compute it using
        % the eigenvalues from a generalized eigenvalue problem
        lam = eig(H,E);
        mask = (lam<0) & (lam>-100*eps(max(abs(lam))));
        lam(mask) = 0;
        L_df1 = p*q;
        L_df2 = r*t-2*u;
        if isreal(lam) && all(lam>=0) && L_df2>0
            L = prod(1./(1+lam));
        else
            L = NaN;
            L_df2 = max(0,L_df2);
        end
        L1 = L^(1/t);
        L_F = ((1-L1) / L1) * (r*t-2*u)/(p*q);
        L_rsq = 1-L1;
        
        Value = L;
        F = L_F;
        RSquare = L_rsq;
        df1 = L_df1;
        df2 = L_df2;
    case 'Pillai'
            error('code that')
    case 'Hotelling'
        lam = eig(H,E);
        mask = (lam<0) & (lam>-100*eps(max(abs(lam))));
        lam(mask) = 0;
        
        % ~~~ Hotelling-Lawley trace = U
        % U = trace( H * E^-1 ) but it can also be written as the sum of the
        % eigenvalues that we already obtained above
        if isreal(lam) && all(lam>=0)
            U = sum(lam);
        else
            U = NaN;
            n(n<0) = NaN;
        end
        b = (p+2*n)*(q+2*n) / (2*(2*n+1)*(n-1));
        c = (2 + (p*q+2)/(b-1))/(2*n);
        if n>0
            U_F = (U/c) * (4+(p*q+2)/(b-1)) / (p*q);
        else
            U_F = U * 2 * (s*n+1) / (s^2 * (2*m+s+1));
        end
        U_rsq = U / (U+s);
        U_df1 = s*(2*m+s+1);
        U_df2 = 2*(s*n+1);
        
        Value = U;
        F = U_F;
        RSquare = U_rsq;
        df1 = U_df1;
        df2 = U_df2;
    case 'Roy'
            error('code that')
    otherwise
        error('x')
end
pValue = fcdf(F, df1, df2, 'upper');
if exist('withinNames','var') && ~isempty(withinNames) && exist('betweenNames','var') && ~isempty(betweenName)
    Within = withinNames;
    Between = betweenNames;
    Statistic = categorical({statLabel}');
    ds = table(Within,Between,Statistic, Value, F, RSquare,df1,df2);
    ds.pValue = pValue;
else
    ds = [];
end

