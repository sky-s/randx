function [R,pdf] = randp(in,N,mFirst,lambda)
%RANDP PERT distributed pseudorandom numbers.
%   R = RANDP([a m b],N) returns an N-by-N matrix containing pseudorandom values
%   drawn from the PERT distribution on the open interval(a,b) with mode
%   m. RANDP([a m b],[M,N]) returns an M-by-N matrix. RANDP([a m b],[M,N,P,...])
%   returns an M-by-N-by-P-by-... array. RANDP([a m b]) returns a scalar. 
% 
%   RANDP([a b],...) uses the midpoint between a and b for m.
%
%   RANDP([a NaN b],...) returns pseudorandom numbers on a uniform distribution
%   between a and b.
% 
%   R = RANDP(m,...), or if a = b, returns scalar R = m.
% 
%   R = RANDP(...,N,mFirst), where mFirst = TRUE, sets R(1) = m.
% 
%   R = RANDP(...,N,mFirst,lambda) reshapes the distribution based on
%   non-negative shape parameter lambda. Default lambda = 4. Higher lambda
%   yields a more "peaky" distribution, lower lambda yields a more uniform
%   distribution (with lambda = 0 being fully uniform). 
% 
%   Note: RANDP takes ~10-50 times longer than randt. Use care with large N.
%
%   Example: Based purely on guesstimates that include a best guess and
%   intuitive upper/lower limits, how much do a dime, nickel, and quarter weigh
%   together?
%       % First, create a function based on RANDP to conveniently define and
%       % generate "uncertain variables."
%       uvar = @(x) randp(x,[1e5,1],true); 
% 
%       dime    = uvar([1 1.5 3.5]); %Pretty light, but not less than a gram..."
%       nickel  = uvar([3 5 6]);
%       quarter = uvar([5 8 10]);
% 
%       total = dime + nickel + quarter;
% 
%       myBestGuess = total(1)
%       % What I would guess if I didn't also capture uncertainty and bracketing
%       % assumptions (true answer is 12.9 grams).
% 
%       histogram(total,'Normalization','pdf');     
% 
%   See also RAND, randn, randt, pctl, betaincinv,
%     invprctile - http://www.mathworks.com/matlabcentral/fileexchange/41131,
%     randraw    - http://www.mathworks.com/matlabcentral/fileexchange/7309.
% 
%   RANDP([a m b],[M,N,...],mFirst,lambda)

%   References: 
%     NASA report: "Analytic Method for Probabilistic Cost and Schedule Risk
%     Analysis" http://goo.gl/dk5d26
%     www.riskamp.com/beta-pert

%   Copyright 2017 Sky Sartorius
%   Contact: www.mathworks.com/matlabcentral/fileexchange/authors/101715

if nargin < 4
    lambda = 4;
end
validateattributes(lambda,{'numeric'},{'scalar','nonnegative','finite'},...
    'randp','shape parameter, lambda,');

if nargin < 3
    mFirst = 0;
end
if nargin < 2
    N = 1;
end

validateattributes(in,{'numeric','DimVar'},{'vector'},'','[a m b] vector')

a = in(1);
b = in(end);

switch numel(in)
    case 3
        m = in(2);
    case {1,2}
        m = (a + b)/2;
    otherwise
        error('Input vector must have 1, 2 ([a b]), or 3 ([a m b]) elements.')
end

if (m>b && m>a) || (m<b && m<a)
    error('Ensure m lies between a and b.')
end

if a == b
    R = m;
    if nargout > 1
        pdf = @(x) (x==a)*1e1000;
    end
    return
end

u = rand(N);

if isnan(m)
    R = a + u*(b-a);
    m = (a + b)/2;
else
    mu = (a + lambda.*m + b) ./ (lambda + 2);
    
    alph = (mu - a).*(2*m - a - b) ./ ((m - mu).*(b - a));
    
    ind = abs(1 - m./mu) < eps; % i.e. m==mu, practically.
    alph(ind) = lambda./2 + 1;
    
    bet = alph.*(b - mu) ./ (mu - a);
  
    R = a + (b - a).*betaincinv(u,alph,bet);
end

if mFirst
    R(1) = m;
end

if nargout > 1
    % Return pdf as function handle.
    pdf = @(x) ( (x-a).^(alph-1)  .*  (b-x).^(bet-1) )   ./...
               ( (b-a).^(alph+bet-1)  .*  beta(alph,bet) );
end

