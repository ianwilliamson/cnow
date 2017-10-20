function A = func_Anorm_dfg(A0,z,freq1,freq2,d)
P10 = A0(1)^2;
P20 = A0(2)^2;
P30 = A0(3)^2;

gamma  = 1j*sqrt((P10/freq1)./(P30/(freq1+freq2)))*ones(1,length(z));
kappa  = d*8.85418782e-12*(freq1+freq2)/2;
GammaL = 1j*kappa*sqrt(freq1/freq2*P30)*z;

P2L = -P10*freq2/freq1.*ellipji(GammaL,gamma).^2;
P1L = P10+freq1/freq2*P2L;
P3L = P30-(freq1+freq2)/freq2*P2L;

A = [sqrt(P1L); sqrt(P2L); sqrt(P3L)];
end