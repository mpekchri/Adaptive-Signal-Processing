% Noise cancellation solving Wiener-Hopf's equation
% author: Christophoros Bekos (8311)

clear all
close all

N = 1000; % number of time steps
M = 10;   % filter's length

%% signal A
%we assume that A is normally distributed
A = sind(1:N);

%% signal x
f = pi/6;
x = zeros(N,1);
for i=1:N
    x(i) = A(i)*sin((pi/8)*N+f);
end

%% white noise N
N_mu = 0;
N_sigma = 0.32;
n = normrnd(N_mu,N_sigma,N,1);

%% input signal u
u = zeros(N,1);
u(1) = n(1);
u(2) = n(2);
for i=3:N
    u(i) = 0.25*u(i-1) - 0.12*u(i-2) + n(i);
end

%% desired signal d
d = zeros(N,1);
for i=1:N
  d(i) = x(i) + n(i);
end

%% some figures
figure
plot([x d])
legend({'x(n)', 'd(n)'})

figure
plot([d u])
legend({'d(n)', 'u(n)'})

%% fir filter
a = xcorr(u,u,M-1,'unbiased');
b = a(M:(2*M-1));
R = toeplitz(b);
a = xcorr(d,u,M,'unbiased');
b = a(M:(2*M));
p = toeplitz(b);
p = p(1,2:end)';
wo = R\p;

% filter apply
y = zeros(N-M,1);
for i = M:N
     y(i) = u(i:-1:i-M+1)' * wo;
end
e = d - y;

figure
plot([y u])
legend({'y(n)', 'u(n)'})

figure
plot([x e])
legend({'x(n)', 'x_estimation(n)'})

%% parameter error
figure
plot((x-e).^2);
title('estimation error')





