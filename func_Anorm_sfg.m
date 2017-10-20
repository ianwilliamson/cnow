function A = func_Anorm_sfg(A0,z,freqS,freqP,d)
P10 = A0(1)^2;
P20 = A0(2)^2;

gamma  = sqrt((P10/freqS)./(P20/freqP));
kappa  = d*8.85418782e-12*(freqS+freqP)/2;
GammaL = kappa*sqrt(freqS/(freqS+freqP)*P20)*z;

P3L = P10*(freqS+freqP)/freqS.*ellipj(GammaL,gamma).^2;
P1L = P10-freqS/(freqS+freqP)*P3L;
P2L = P20-freqP/(freqS+freqP)*P3L;

A = [sqrt(P1L); sqrt(P2L); sqrt(P3L)];
end