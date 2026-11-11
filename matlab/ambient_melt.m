function mdot = ambient_melt(T,S,z,v,p)
% AMBIENT_MELT  Melt rate (m/d) from the three-equation formulation, no plume.
%   T = conservative temperature (degC), S = absolute salinity (g/kg),
%   z = depth (m, NEGATIVE down), v = background velocity (m/s), p = parameter struct.
%   T, S and z are column vectors of equal length.
lambda1=p.lambda1; lambda2=p.lambda2; lambda3=p.lambda3;
GT=p.GT; GS=p.GS; cw=p.cw; ci=p.ci; L=p.L; k=p.k; Tice=p.Tice;

a = lambda1*(GT*cw-GS*ci);
b = GS*ci*(lambda1*S-lambda2-lambda3*z+Tice-(L/ci)) - GT*cw*(T-lambda2-lambda3*z);
c = GS*S.*(ci*(lambda2+lambda3*z-Tice)+L);
Sb = (1./(2*a)).*(-b-sqrt(b.^2-4*a*c));
mdot = 86400*GS*sqrt(k)*v.*(S-Sb)./Sb;     % m/d
end