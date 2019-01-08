% NOTE: Based on Matlab code in Olivier Ledoit and Michael Wolf. Analytical
%       Nonlinear Shrinkage of Large-Dimensional Covariance Matrices. (Nov 2018)
%       http://www.econ.uzh.ch/static/wp/econwp264.pdf
% NOTE: this assumes that X has centered columns
function sigmatilde=analytical_shrinkage(X, lambda, u)
    % extract sample eigenvalues sorted in ascending order and eigenvectors
    [n, p] = size(X);
    meanx=mean(X);
    X=X-meanx(ones(n,1),:);
    % important: sample size n must be >= 12
    sample = (X' * X) ./ n;
    #[u, lambda] = eig(sample, 'vector');
    [lambda, isort] = sort(lambda);
    u = u(:, isort);
    % compute analytical nonlinear shrinkage kernel formula
    lambda = lambda(max(1,p-n+1):p);
    L = repmat(lambda,[1 min(p,n)]);
    h = n^(-1/3);
    H = h*L';
    x = (L-L')./H;
    ftilde = (3/4/sqrt(5))*mean(max(1-x.^2./5,0)./H,2);
    Hftemp = (-3/10/pi)*x+(3/4/sqrt(5)/pi)*(1-x.^2./5).*log(abs((sqrt(5)-x)./(sqrt(5)+x)));
    Hftemp(abs(x)==sqrt(5)) = (-3/10/pi)*x(abs(x)==sqrt(5));
    Hftilde = mean(Hftemp./H, 2);
    if p<=n
      dtilde = lambda./((pi*(p/n)*lambda.*ftilde).^2+(1-(p/n)-pi*(p/n)*lambda.*Hftilde).^2);
    else
      Hftilde0 = (1/pi)*(3/10/h^2+3/4/sqrt(5)/h*(1-1/5/h^2)*log((1+sqrt(5)*h)/(1-sqrt(5)*h)))*mean(1./lambda);
      dtilde0 = 1/(pi*(p-n)/n*Hftilde0);
      dtilde1 = lambda./(pi^2*lambda.^2.*(ftilde.^2+Hftilde.^2));
      dtilde = [dtilde0*ones(p-n,1);dtilde1];
    end
    sigmatilde=u*diag(dtilde)*u';
end

C_20  = dlmread("tmpmat/C_20.csv");
C_40  = dlmread("tmpmat/C_40.csv");
C_200 = dlmread("tmpmat/C_200.csv");
C_400 = dlmread("tmpmat/C_400.csv");

X_40x20   = dlmread('tmpmat/X_40x20.csv');
X_20x40   = dlmread('tmpmat/X_20x40.csv');
X_400x200 = dlmread('tmpmat/X_400x200.csv');
X_200x400 = dlmread('tmpmat/X_200x400.csv');

l_40x20 = dlmread('tmpmat/S_lambda_40x20.csv');
l_20x40 = dlmread('tmpmat/S_lambda_20x40.csv');
l_200x400 = dlmread('tmpmat/S_lambda_200x400.csv');
l_400x200 = dlmread('tmpmat/S_lambda_400x200.csv');

u_40x20 = dlmread('tmpmat/S_u_40x20.csv');
u_20x40 = dlmread('tmpmat/S_u_20x40.csv');
u_200x400 = dlmread('tmpmat/S_u_200x400.csv');
u_400x200 = dlmread('tmpmat/S_u_400x200.csv');


C = {C_20, C_40, C_200, C_400};
X = {X_40x20, X_20x40, X_400x200, X_200x400};
L = {l_40x20, l_20x40, l_400x200, l_200x400};
U = {u_40x20, u_20x40, u_400x200, u_200x400};

times = zeros(4, 1)
res = zeros(4, 1)

for (i = [1, 3])
  Xi = X{1, i};
  Ci = C{1, i};
  Li = L{1, i};
  Ui = U{1, i};
  tic();
  C_est = analytical_shrinkage(Xi, Li, Ui);
  times(i) = toc();
  res(i) = sqrt(sum(sum((Ci - C_est).^2)));
end
