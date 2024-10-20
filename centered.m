%% How to shift the spectrum to have zero-frequency at the center
%%
%%    fft(X) : Yk = somme(i=[0,n[, Xi * exp(-2jpi * k * i / n)
%%
%%     let Zi = Xi * exp(2jpi * i * x/n)
%%
%%    fft(Z) : somme(i=[0,n[, Xi * exp(-2jpi * k * (i-x) / n)
%%    => Using x = n/2 to have fft(Z) = fftshift(fft(X))
%%             Zi = Xi * (-1 if odd(i))


% Number of points
n = 128;
t = transpose(0:(n-1));

% Signal
s = rand(n,1) + 1i*rand(n,1);

% FFT with continuous-F at center
r = fftshift(fft(s));

% FFT of signal multiplied by appropriate frequency to center continuous-F
d = exp(-2*pi*1i* t/2);   % i.e. :           Here error-max=1e-13
d = ones(n,1); d(2:2:end) = -d(2:2:end);   % Then error-max=0
f = fft(s .* d);

plot(t, r, t, f+10);

fprintf("Error-max = %e\n", max(abs(f-r)));
