% function I=supermarcumQ(a,b,errabs)
%
% Efficient computation of Marcum-Q function
%
% supermarcumq(a,b) --> evaluates Marcum-Q function with 1e-10 error
% tolerance
%
% a,b non-negative real values
%
% errabs: absolute error (OPTIONAL) (DEFAULT 1e-10)


function I=supermarcumq(a,b,errabs)

if (nargin==2),
errabs=1e-10;
end


if (b>a),

z=a./b;

I=1/pi*integral(@integrando1,0,pi);

elseif (b<a)

z=b./a;

I=1+1/pi*integral(@integrando2,0,pi);

else 

I=(1+besseli(0,a.^2,1))/2;

end

function y = integrando1(t)

y = (1+z.*cos(t))./(1+2.*z.*cos(t)+z.^2) .*exp( -b.^2/2.*(1+2*z.*cos(t)+z.^2) );

end

function y = integrando2(t)

y = (z.^2+z.*cos(t))./(1+2.*z.*cos(t)+z.^2) .*exp( -a.^2/2.*(1+2*z.*cos(t)+z.^2) );

end

end