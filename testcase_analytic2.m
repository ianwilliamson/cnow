fact0 = 1999.99993625136;

d     = 10e-4;
L     = 0.075e-2;
z     = linspace(0,L,99);
Delta = 1e9;
freqS = 190e12;
freqP = 210e12;
omega = [freqS,freqP,freqS+freqP];

k = @(p,zi) omega(p)/3e8;

A0 = [0,10,1];

%%
cnow = Coupled_Nonlinear_Optical_Waves(d,z,A0,k,omega);
cnow.solve()
A_newton = abs(cnow.A);
A_analytic = func_Anorm_dfg(A0,z,freqS,freqP,d*fact0);

figure(); hold on; box on; grid on;
plot(z,A_analytic,'ko-');
plot(z,A_newton,'ro-');
set(gca(),'yscale','log');
ylabel('|A_n|');
xlim([min(z),max(z)]);

error_norm = norm(A_newton-A_analytic)

%% Fitting
func_fit = @(factor) norm(A_newton-func_Anorm_sfg(A0,z,freqS,freqP,d*factor));
fact0 = fminbnd(func_fit,1000,2000)