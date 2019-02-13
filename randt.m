function R = randt(in,N,mFirst)
%RANDT Triangularly distributed pseudorandom numbers.
%   R = RANDT([a m b],N) returns an N-by-N matrix containing pseudorandom values
%   drawn from the triangular distribution on the open interval(a,b) with mode
%   m. RANDT([a m b],[M,N]) returns an M-by-N matrix. RANDT([a m b],[M,N,P,...])
%   returns an M-by-N-by-P-by-... array. RANDT([a m b]) returns a scalar. 
% 
%   RANDT([a b],...) uses the midpoint between a and b for m.
%
%   RANDT([a NaN b],...) returns pseudorandom numbers on a uniform distribution
%   between a and b.
% 
%   R = RANDT(m,...), or if a = b, returns scalar R = m.
% 
%   R = RANDT(...,N,mFirst), where mFirst = TRUE, sets R(1) = m.
%
%   Example: Based purely on guesstimates that include a best guess and
%   intuitive upper/lower limits, how much do a dime, nickel, and quarter weigh
%   together?
%       % First, create a function based on RANDT to conveniently define and
%       % generate "uncertain variables."
%       uvar = @(x) randt(x,[1e5,1],true); 
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
%   See also RAND, randn, pctl, 
%     invprctile - http://www.mathworks.com/matlabcentral/fileexchange/41131,
%     randraw    - http://www.mathworks.com/matlabcentral/fileexchange/7309.
% 
%   RANDT([a m b],[M,N,...],mFirst)

%   Copyright 2017 Sky Sartorius
%   Contact: www.mathworks.com/matlabcentral/fileexchange/authors/101715

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
    return
end

u = rand(N);

if isnan(m)
    R = a + u*(b-a);
    m = (a + b)/2;
else
    t = (m-a) / (b-a);
    R = a + (b-a)*((u <= t).*sqrt(t*u) + (u > t).*(1-sqrt((1-t)*(1-u))));
end

if mFirst
    R(1) = m;
end