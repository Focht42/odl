function ans = bessj0(x)
%% Bessel function of order zero, which is much faster than matlab's besselj(0,x)
%
% matlab translation of a C code in
% Cephes Math Library Release 2.1:  January, 1989
% Copyright 1984, 1987, 1989 by Stephen L. Moshier
%
% The domain is divided into the intervals [0, 5] and
% (5, infinity). In the first interval the following rational
% approximation is used:
%
%
%        2         2
% (w - r  ) (w - r  ) P (w) / Q (w)
%       1         2    3       8
%
%            2
% where w = x  and the two r's are zeros of the function.
%
% In the second interval, the Hankel asymptotic expansion
% is employed with two rational functions of degree 6/6
% and 7/7.
%
%
% ACCURACY:
%
%                      Absolute error:
% arithmetic   domain     # trials      peak         rms
%    IEEE      0, 30       60000       4.2e-16     1.1e-16

SQ2OPI =  7.9788456080286535587989E-1;  % sqrt( 2/pi )
PIO4   =  7.85398163397448309616E-1;    % pi/4

PP = [7.96936729297347051624E-4;...
    8.28352392107440799803E-2;...
    1.23953371646414299388E0;...
    5.44725003058768775090E0;...
    8.74716500199817011941E0;...
    5.30324038235394892183E0;...
    9.99999999999999997821E-1];
PQ = [9.24408810558863637013E-4;...
    8.56288474354474431428E-2;...
    1.25352743901058953537E0;...
    5.47097740330417105182E0;...
    8.76190883237069594232E0;...
    5.30605288235394617618E0;...
    1.00000000000000000218E0];

RP = [-4.79443220978201773821E9;...
    1.95617491946556577543E12;...
    -2.49248344360967716204E14;...
    9.70862251047306323952E15];

RQ = [...% 1.00000000000000000000E1
    4.99563147152651017219E2;...
    1.73785401676374683123E5;...
    4.84409658339962045305E7;...
    1.11855537045356834862E10;...
    2.11277520115489217587E12;...
    3.10518229857422583814E14;...
    3.18121955943204943306E16;...
    1.71086294081043136091E18];


QP = [-1.13663838898469149931E-2;...
    -1.28252718670509318512E0;...
    -1.95539544257735972385E1;...
    -9.32060152123768231369E1;....
    -1.77681167980488050595E2;...
    -1.47077505154951170175E2;...
    -5.14105326766599330220E1;...
    -6.05014350600728481186E0];
QQ = [...%1.00000000000000000000E0;...
    6.43178256118178023184E1;...
    8.56430025976980587198E2;...
    3.88240183605401609683E3;...
    7.24046774195652478189E3;...
    5.93072701187316984827E3;...
    2.06209331660327847417E3;...
    2.42005740240291393179E2];

DR1 = 5.78318596294678452118E0;
DR2 = 3.04712623436620863991E1;

ans = zeros(size(x));
ind = find(x<0);
if ~isempty(ind)
x(ind) = -x(ind)
end

ind1 = find(x<1e-5);
ind2 = find(1e-5<=x & x<=5.0);
ind3 = find(x>5.0);

if ~isempty(ind1)
    z= x(ind1);
    ans(ind1) = 1-z.*z/4;
end
if ~isempty(ind2)
    z = x(ind2);
    z = z.*z;
    p = (z - DR1).*(z - DR2);
    ans(ind2) = p .* polevl(z, RP)./p1evl(z, RQ);
end
if ~isempty(ind3)
    one = ones(size(ind3));
    xx = x(ind3);
    w = 5.0*one./xx;
    q = 25.0 *one./(xx.*xx);
    p = polevl(q, PP) ./ polevl(q, PQ);
    q = polevl(q, QP) ./ p1evl(q,  QQ);
    xn = xx - PIO4;
    p = p.* cos(xn) - w.*q.* sin(xn);
    ans(ind3) = SQ2OPI*p./ sqrt(xx);
end
end

function ans=polevl(x,coef)
%% evaluates polynomial coeff(1)x^N+..+coef(N+1)
ans = coef(1) * ones(size(x));
for j=2:length(coef)
    ans = ans .* x  + coef(j);
end
end

function ans= p1evl(x, coef)
%% evaluates polynomial x^N + coef(1)*x^(N-1) + ..+ coef(N)
ans = x + coef(1);
for j=2:length(coef)
    ans = ans.* x  + coef(j);
end
end


