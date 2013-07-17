function [rho,al0,lambda,p0] = density_wright_eos(T, S, P)
% density_wright_eos  Density of seawater using Wright EOS
%
% density_wright_eos(T,S,P)
%
% Returns the in situ density in kg/m^3.
%
% T - potential temperature relative to the surface in C.
% S - salinity in PSU.
% P - pressure in Pa.
%
% Reference: Wright, Daniel, 1997.
% http://dx.doi.org/10.1175/1520-0426(1997)014%3C0735:AEOSFU%3E2.0.CO;2
%
% Coded for Matlab by W. Anderson 9/07
% Updated for gtoosl by AJA, Fall 2011

if nargin ~= 3
   error('density_wright_eos.m: Must pass 3 parameters ')
end %if

if min(S(:))<0
   error('Salinity < 0! You might have the order of arguments wrong...')
end

% Extended formula fit over reduced range (Table 1,column 4)
%a0 = 7.057924e-4; a1 = 3.480336e-7; a2 = -1.112733e-7;
%b0 = 5.790749e8;  b1 = 3.516535e6;  b2 = -4.002714e4;
%b3 = 2.084372e2;  b4 = 5.944068e5;  b5 = -9.643486e3;
%c0 = 1.704853e5;  c1 = 7.904722e2;  c2 = -7.984422;
%c3 = 5.140652e-2; c4 = -2.302158e2; c5 = -3.079464;

% Extended formula fit over full range (Table 1,column 2)
a0 = 7.133718e-4; a1 = 2.724670e-7; a2 = -1.646582e-7;
b0 = 5.613770e8;  b1 = 3.600337e6;  b2 = -3.727194e4;
b3 = 1.660557e2;  b4 = 6.844158e5;  b5 = -8.389457e3;
c0 = 1.609893e5;  c1 = 8.427815e2;  c2 = -6.931554;
c3 = 3.869318e-2; c4 = -1.664201e2; c5 = -2.765195;

al0 = a0 + (a1*T +a2*S);
p0 = (b0 + b4*S) + T .* (b1 + T.*(b2 + b3*T) + b5*S);
lambda = (c0 +c4*S) + T .* (c1 + T.*(c2 + c3*T) + c5*S);
rho = (P + p0) ./ (lambda + al0.*(P + p0));
