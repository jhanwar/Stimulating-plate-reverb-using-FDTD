clear variables;
%% Params
SR = 44100;
k = 1 / SR;
simTime = 0.01;
T = SR * simTime;
showAnimation = 1;
kappa = 1;
gamma = 50;
s0 = 700;

% [F, SR] = audioread('Somewhere In My Car.mp3');
% T = length(F);
% k = 1 / SR;
% F = SR^2 * F;
% Comment out below and uncomment above for audio input
F = zeros(1, T);
F(1) = 5e8;

N = floor(sqrt((-gamma^2 / kappa^2 + sqrt(gamma^4 / kappa^4 + 16 / (k^2 * kappa^2))) / 16));
h = 1 / N;
m2 = kappa^2 * k^2 / h^4;
l2 = gamma^2 * k^2 / h^2;

%% Setup
up = zeros(N, N);
u = zeros(N, N);
un = zeros(N, N);

% u = sin(pi * (0 : N - 1)/(N - 1)).^2' * sin(pi * (0 : N - 1)/(N - 1)).^2;
% u(floor(N/2) - floor(N/6) : floor(N/2) + floor(N/6), floor(N/2) - floor(N/6) : floor(N/2) + floor(N/6)) = 1;
% u = 1 * sin(pi * (0 : N - 1)/(N - 1))' * sin(pi * (0 : N - 1)/(N - 1));
% u = rand(N, N);
cx = floor(N / 2);
cy = floor(N / 2);
hw = 4;
% u(cy - hw : cy + hw, cx - hw : cx + hw) = sin(pi * (0 : 2 * hw) / (2 * hw)).^2' * sin(pi * (0 : 2 * hw) / (2 * hw)).^2;
% up = u;
e = zeros(N, N);
[xfac, yfac] = meshgrid(sin(pi * (0 : 2 * hw) / (2 * hw)).^2, sin(pi * (0 : 2 * hw) / (2 * hw)).^2);
e(cy - hw : cy + hw, cx - hw : cx + hw) = xfac .* yfac;
% e(cx, cy) = 1;
y = zeros(1, T);
K = zeros(1, T);
P = zeros(1, T);

%% FDTD Scheme
if showAnimation
    frames(T) = struct('cdata', [], 'colormap', []); zmax = 1.5;
    figure;
    s = surf((0 : N - 1) / (N - 1), (0 : N - 1) / (N - 1), u);
    xlabel('X');
    ylabel('Y');
    zlabel('Displacement');
    title('Vibration of stiff plate under tension');
    set(gca, 'xtick', 0 : 0.2 : 1, 'ytick', 0 : 0.2 : 1);
    set(s, 'edgecolor', 'none', 'FaceColor', 'white')
    material shiny
    camlight
    zlim([-zmax zmax]);
end
for n = 1 : T
    %     for l = 3 : N - 2
    %         for m = 3 : N - 2
    %             un(l, m) = (2 - 20 * m2) * u(l, m) + 8 * m2 * (u(l, m + 1) + u(l, m - 1) + u(l + 1, m) + u(l - 1, m)) ...
    %                 - 2 * m2 * (u(l + 1, m + 1) + u(l + 1, m - 1) + u(l - 1, m + 1) + u(l - 1, m - 1)) ...
    %                 - m2 * (u(l, m + 2) + u(l, m - 2) + u(l + 2, m) + u(l - 2, m)) - up(l, m);
    %         end
    %     end
    un(3 : N - 2, 3 : N - 2) = ((2 - 20 * m2 - 4 * l2) * u(3 : N - 2, 3 : N - 2) ...
        + (8 * m2 + l2) * (u(3 : N - 2, 4 : N - 1) + u(3 : N - 2, 2 : N - 3) + u(4 : N -1, 3 : N - 2) + u(2 : N - 3, 3 : N - 2)) ...
        - 2 * m2 * (u(4 : N - 1, 4 : N - 1) + u(4 : N - 1, 2 : N - 3) + u(2 : N - 3, 4 : N - 1) + u(2 : N - 3, 2 : N - 3)) ...
        - m2 * (u(3 : N - 2, 5 : N) + u(3 : N - 2, 1 : N - 4) + u(5 : N, 3 : N - 2) + u(1 : N - 4, 3 : N - 2)) - (1 - s0 * k) * up(3 : N - 2, 3 : N - 2)) / (1 + s0 * k);
    
    % Vertical Edge
    %     for l = 3 : N - 2
    %         un(l, 2) = (2 - 20 * m2) * u(l, 2) + 8 * m2 * (u(l, 3) + u(l + 1, 2) + u(l - 1, 2)) ...
    %             - 2 * m2 * (u(l + 1, 3) + u(l - 1, 3)) ...
    %             - m2 * (u(l, 4) + u(l + 2, 2) + u(l - 2, 2)) - up(l, 2);
    %         un(l, N - 1) = (2 - 20 * m2) * u(l, N - 1) + 8 * m2 * (u(l, N - 2) + u(l + 1, N - 1) + u(l - 1, N - 1)) ...
    %             - 2 * m2 * (u(l + 1, N - 2) + u(l - 1, N - 2)) ...
    %             - m2 * (u(l, N - 3) + u(l + 2, N - 1) + u(l - 2, N - 1)) - up(l, N - 1);
    %     end
    un(3 : N - 2, 2) = ((2 - 20 * m2 - 4 * l2) * u(3 : N - 2, 2) ...
        + (8 * m2 + l2) * (u(3 : N - 2, 3) + u (4 : N - 1, 2) + u(2 : N - 3, 2)) ...
        - 2 * m2 * (u(4 : N - 1, 3) + u(2 : N - 3, 3)) ...
        - m2 * (u(3 : N - 2, 4) + u(5 : N, 2) + u(1 : N - 4, 2)) - (1 - s0 * k) * up(3 : N - 2, 2)) / (1 + s0 * k);
    un(3 : N - 2, N - 1) = ((2 - 20 * m2 - 4 * l2) * u(3 : N - 2, N - 1) ...
        + (8 * m2 + l2) * (u(3 : N - 2, N - 2) + u (4 : N - 1, N - 1) + u(2 : N - 3, N - 1)) ...
        - 2 * m2 * (u(4 : N - 1, N - 2) + u(2 : N - 3, N - 2)) ...
        - m2 * (u(3 : N - 2, N - 3) + u(5 : N, N - 1) + u(1 : N - 4, N - 1)) - (1 - s0 * k) * up(3 : N - 2, N - 1)) / (1 + s0 * k);
    
    % Horizontal Edge
    %     for m = 3 : N - 2
    %         un(2, m) = (2 - 20 * m2) * u(2, m) + 8 * m2 * (u(3, m) + u(2, m + 1) + u(2, m - 1)) ...
    %             - 2 * m2 * (u(3, m + 1) + u(3, m - 1)) ...
    %             - m2 * (u(4, m) + u(2, m + 2) + u(2, m - 2)) - up(2, m);
    %         un(N - 1, m) = (2 - 20 * m2) * u(N - 1, m) + 8 * m2 * (u(N - 2, m) + u(N - 1, m + 1) + u(N - 1, m - 1)) ...
    %             - 2 * m2 * (u(N - 2, m + 1) + u(N - 2, m - 1)) ...
    %             - m2 * (u(N - 3, m) + u(N - 1, m + 2) + u(N - 1, m - 2)) - up(N - 1, m);
    %     end
    un(2, 3 : N - 2) = ((2 - 20 * m2 - 4 * l2) * u(2, 3 : N - 2) ...
        + (8 * m2 + l2) * (u(3, 3 : N - 2) + u(2, 4 : N - 1) + u(2, 2 : N - 3)) ...
        - 2 * m2 * (u(3, 4 : N - 1) + u(3, 2 : N - 3)) ...
        - m2 * (u(4, 3 : N - 2) + u(2, 5 : N) + u(2, 1 : N - 4)) - (1 - s0 * k) * up(2, 3 : N - 2)) / (1 + s0 * k);
    un(N - 1, 3 : N - 2) = ((2 - 20 * m2 - 4 * l2) * u(N - 1, 3 : N - 2) ...
        + (8 * m2 + l2) * (u(N - 2, 3 : N - 2) + u(N - 1, 4 : N - 1) + u(N - 1, 2 : N - 3)) ...
        - 2 * m2 * (u(N - 2, 4 : N - 1) + u(N - 2, 2 : N - 3)) ...
        - m2 * (u(N - 3, 3 : N - 2) + u(N - 1, 5 : N) + u(N - 1, 1 : N - 4)) - (1 - s0 * k) * up(N - 1, 3 : N - 2)) / (1 + s0 * k);
    
    % Corners
    un(2, 2) = ((2 - 20 * m2 - 4 * l2) * u(2, 2) + (8 * m2 + l2) * (u(2, 3) + u(3, 2)) ...
        - 2 * m2 * (u(3, 3)) ...
        - m2 * (u(2, 4) + u(4, 2)) - (1 - s0 * k) * up(2, 2)) / (1 + s0 * k);
    un(2, N - 1) = ((2 - 20 * m2 - 4 * l2) * u(2, N - 1) + (8 * m2 + l2) * (u(2, N - 2) + u(3, N - 1)) ...
        - 2 * m2 * (u(3, N - 2)) ...
        - m2 * (u(2, N - 3) + u(4, N - 1)) - (1 - s0 * k) * up(2, N - 1)) / (1 + s0 * k);
    un(N - 1, 2) = ((2 - 20 * m2 - 4 * l2) * u(N - 1, 2) + (8 * m2 + l2) * (u(N - 2, 2) + u(N - 1, 3)) ...
        - 2 * m2 * (u(N - 2, 3)) ...
        - m2 * (u(N - 3, 2) + u(N - 1, 4)) - (1 - s0 * k) * up(N - 1, 2)) / (1 + s0 * k);
    un(N - 1, N - 1) = ((2 - 20 * m2 - 4 * l2) * u(N - 1, N - 1) + (8 * m2 + l2) * (u(N - 2, N - 1) + u(N - 1, N - 2)) ...
        - 2 * m2 * (u(N - 2, N - 2)) ...
        - m2 * (u(N - 3, N - 1) + u(N - 1, N - 3)) - (1 - s0 * k) * up(N - 1, N - 1)) / (1 + s0 * k);
    
    % Forcing term
    un = un + k^2 * F(n) * e / (1 + s0 * k);
    
    KE = 0.5 * sum(sum((u - up).^2 / k^2)) * h^2;
    PE = kappa^2 * 0.5 * sum(sum((u(3 : N, 2 : N - 1) + u(1 : N - 2, 2 : N - 1) + u(2 : N - 1, 3 : N) + u(2 : N - 1, 1 : N - 2) - 4 * u(2 : N - 1, 2 : N - 1))/(h^2) ...
        .* (up(3 : N, 2 : N - 1) + up(1 : N - 2, 2 : N - 1) + up(2 : N - 1, 3 : N) + up(2 : N - 1, 1 : N - 2) - 4 * up(2 : N - 1, 2 : N - 1))/(h^2))) * h^2;
    up = u;
    u = un;
    y(n) = u(floor(3 * N / 4), floor(3 * N / 4));
    K(n) = KE;
    P(n) = PE;
    
    if showAnimation
        s.ZData = u;
        frames(n) = getframe;
    end
end

subplot(311);
plot((0:T-1) * k, y);
xlabel('time');
ylabel('displacement');
title('displacement at pickup point');
grid on;

subplot(312);
plot((0:T-1) * k, K);
hold on;
plot((0:T-1) * k, P);
plot((0:T-1) * k, K + P);
legend({'Kinetic', 'Potential', 'Total'});
xlabel('time');
ylabel('numerical energy');
title('numerical energies');
grid on;

subplot(313);
plot((0:T-1) * k, ((K + P) - (K(1) + P(1)))/(K(1) + P(1)));
xlabel('time');
ylabel('relative error');
title('relative error in total numerical energy');
grid on;