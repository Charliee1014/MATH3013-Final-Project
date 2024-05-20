% Define parameters
clc;
clear;
L = 1.7;          % Length of the tub
W = 1;            % Width of the tub
H = 0.7;          % Height of the tub
h_t = 40;         % Heat transfer coefficient (top)
h_b = 16;         % Heat transfer coefficient (bottom)
h_w = 32;         % Heat transfer coefficient (sides)
hS = (h_t + h_b) * L * W + h_w * (2 * H * W + 2 * H * L); % Total heat transfer coefficient
c = 4178;         % Specific heat capacity of water
m_tub = 1190;     % Mass of water in the tub
T_inf = 293.15;   % Ambient temperature
T0 = 310.15;      % Initial water temperature
T_in = 318.15;    % Input water temperature

% Calculate lambda
lambda = hS / (c * m_tub);

% Define the function f1 for temperature variation over time
f1 = @(t) T_inf + (T0 - T_inf) * exp(-lambda * t);

% Time range from 0 to 1800 seconds
t1 = 0:1:1800;

% Calculate the temperature at different time points
y1 = f1(t1) - 273.15;

% Plot the graph
figure;
plot(t1, y1);
xlabel('Time (s)');
ylabel('Temperature (°C)');
title('Temperature of Water in Tub Over Time');
grid on;

%%
% Temperature change before and after turning on hot water
m_dot = 0.4;
alpha = (-hS - m_dot * c) / (c * m_tub);
beta = (hS * T_inf + m_dot * c * T_in) / (c * m_tub);
Temp_switch = f1(1800);  
C1 = (Temp_switch + beta / alpha);  

% Update the definition of function f2
f2_new = @(t) -beta / alpha + C1 * exp(alpha * (t - 1800));

% Time range
t2 = 1800:1:3600;

% Calculate the temperature
y2 = f2_new(t2) - 273.15; 

% Plot the graph
figure;
hold on;
plot(t1, y1, 'b');  
plot(t2, y2, 'r');  
plot(1800, Temp_switch - 273.15, 'ko', 'MarkerSize', 10);  
text(1800, Temp_switch - 273.15, ' Hot water turned on', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
xlabel('Time (s)');
ylabel('Temperature (°C)');
title('Temperature of Water in Tub Over Time');
grid on;
legend('Before turning on the hot water', 'After turning on the hot water');
hold off;

%%
% Effect of mass flow rate on steady temperature
m_dot = linspace(0, 0.5, 100);  
k = zeros(size(m_dot));  
for i = 1:length(m_dot)
    alpha = (-hS - m_dot(i) * c) / (c * m_tub);
    beta = (hS * T_inf + m_dot(i) * c * T_in) / (c * m_tub);
    k(i) = -beta / alpha - 273.15;  
end

% Plot the graph
figure;
plot(m_dot, k);
xlabel('Mass flow rate (kg/s)');
ylabel('Steady Water Temperature (°C)');
title('Effect of Different Mass Flow Rate');
grid on;

%%
% Effect of input water temperature on steady temperature
T_in_C = linspace(20, 50, 100);  
T_in_K = T_in_C + 273.15;        
k = zeros(size(T_in_C));  
for i = 1:length(T_in_K)
    alpha = (-hS - m_dot * c) / (c * m_tub);
    beta = (hS * T_inf + m_dot * c * T_in_K(i)) / (c * m_tub);
    k(i) = -beta / alpha - 273.15;  
end

% Plot the graph
figure;
plot(T_in_C, k);
xlabel('Input Water Temperature (°C)');
ylabel('Steady Water Temperature (°C)');
title('Effect of Different Input Water Temperature');
grid on;

%%
% Effect of ambient temperature on steady temperature
T_inf_C = linspace(20, 40, 100);  
T_inf_K = T_inf_C + 273.15;       
k = zeros(size(T_inf_C));  
for i = 1:length(T_inf_K)
    beta = (hS * T_inf_K(i) + m_dot * c * T_in) / (c * m_tub);
    k(i) = -beta / alpha - 273.15; 
end

% Plot the graph of k 
figure;
plot(T_inf_C, k);
xlabel('Ambient Air Temperature (°C)');
ylabel('Steady Water Temperature (°C)');
title('Effect of Different Ambient Air Temperature');
grid on;

%% 
% Flow rate 0.1
m_dot = 0.1;
alpha = (-hS - m_dot*c)/(c*m_tub);
beta = (hS*T_inf + m_dot*c*T_in)/(c*m_tub);
k = -beta/alpha - 273.15;
C = 310.15 + beta/alpha;
f3 = @(t) -beta/alpha + C*exp(alpha*t);
t = 0:1:3600;
y3 = f3(t) - 273.15;
plot(t, y3)
xlabel('Time (s)')
ylabel('Temperature (°C)')
title('Temperature of Water in Tub Over Time')
grid on

% Add text
text(2600, 36.94, '$$\dot{m} = 0.1 \, kg/s$$', 'Interpreter', 'latex', 'FontSize', 15)

%% 
clear;
clc;

% Parameters
L = 1.7; % Tub length, m
W = 1; % Tub width, m
H = 0.7; % Tub height, m
rho = 1000; % Water density, kg/m^3
c = 4178; % Specific heat capacity of water, J/(kg*K)
h_w = 32; % Heat transfer coefficient of the wall, W/(m^2*K)
h_b = 16; % Heat transfer coefficient of the bottom, W/(m^2*K)
h_t = 40; % Heat transfer coefficient of the top, W/(m^2*K)
u_inf = 20; % Ambient temperature, °C
alpha = 0.6 / (rho * c); % Thermal diffusivity, m^2/s
Nx = 300; % Number of spatial steps
Nt = 3600; % Number of time steps
dx = L / Nx; % Spatial step size
dt = 1; % Time step size, s
x = linspace(0, L, Nx); % Spatial grid
t = linspace(0, Nt*dt, Nt); % Time grid, from small to large

% Initialize temperature matrix
U = zeros(Nx, Nt);

% Initial conditions
U(:, 1) = 37; % Initial temperature distribution

% Calculate h_tilde
h_tilde = (2 * h_w / W + h_b / H + h_t / H) / (rho * c);


% Time evolution
for n = 1:Nt-1
    for i = 2:Nx-1
        U(i, n+1) = U(i, n) + alpha * dt / dx^2 * (U(i+1, n) - 2*U(i, n) + U(i-1, n)) - h_tilde * dt * (U(i, n) - u_inf);
    end
    U(1, n+1) = U(2, n+1); % Neumann boundary condition, left boundary
    U(Nx, n+1) = U(Nx-1, n+1); % Neumann boundary condition, right boundary
end

% Plot the graph
figure;
imagesc(x, t, U.'); % Use transposed U and adjusted time vector t
colorbarHandle = colorbar; % Get the handle to the colorbar
ylabel(colorbarHandle, 'Temperature (°C)') % Add label to the colorbar
xlabel('Position (m)');
ylabel('Time (s)');
title('Temperature distribution over time and position');

% Set y-axis ticks and labels
yticks = linspace(1, Nt, 11); % Select 11 tick points
yticks = round(yticks); % Ensure yticks are integers
yticklabels = arrayfun(@(y) sprintf('%.0f', t(y)), yticks, 'UniformOutput', false);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);

% Ensure y-axis direction is correct
set(gca, 'YDir', 'normal');

clc;
clear;
bathtubPDE

function bathtubPDE
    % Parameters
    L = 1.7; % Length of the bathtub
    u_in = 45; % Input water temperature
    u_inf = 20; % Ambient temperature
    W = 1; H = 0.7; % Dimensions of the bathtub
    mdot = 0.4; k = 0.6; c = 4178; % Physical parameters
    h_t = 40; h_b = 16; h_w = 32; % Heat transfer coefficients
    alpha = 1.433e-7; % Thermal diffusivity
    h_hat = mdot * c / (k * W * H);
    h_tilde = 1 / (1000 * c) * (2 * h_w / W + h_b / H + h_t / H);

    m = 0; 
    xmesh = linspace(0, L, 600); % Spatial grid
    tspan = linspace(0, 3600, 3600); % Time grid
    options = odeset('RelTol',1e-4,'AbsTol',1e-7); 

    % Solving the PDE for the first part
    sol = pdepe(m, @pdefun, @icfun, @bcfun, xmesh, tspan, options);

    % Plotting the results for the first part
    figure;
    imagesc(xmesh, tspan, sol);
    title('Temperature Distribution Over Time and Position')
    xlabel('Position (m)')
    ylabel('Time (s)')
    colorbarHandle = colorbar;
    ylabel(colorbarHandle, 'Temperature (°C)') 
    axis xy; 

    % Parameters for the second part
    L = 0.1; % Length of the bathtub
    xmesh = linspace(0, L, 600); 
    tspan = linspace(0, 2000, 2000);
    alpha_values = [1e-7, 1.433e-7, 2e-7, 2.5e-7, 3e-7];
    colors = ['r', 'g', 'b', 'c', 'm'];

    figure;
    hold on;

    for i = 1:length(alpha_values)
        alpha = alpha_values(i);
        % Solving the PDE for different alpha values
        sol = pdepe(m, @pdefun, @icfun, @bcfun, xmesh, tspan, options);

        % Plotting the results at t = 2000s
        plot(xmesh, sol(end, :), colors(i), 'DisplayName', ['\alpha = ', num2str(alpha)]);
    end

    title('Temperature Distribution at t = 2000s for Different \alpha Values')
    xlabel('Position (m)')
    ylabel('Temperature (°C)')
    legend('show')
    grid on
    hold off;

    function [c, f, s] = pdefun(x, t, u, DuDx)
        c = 1;
        f = alpha * DuDx;
        s = h_tilde * (u_inf - u);
    end

    function u0 = icfun(x)
        u0 = 37; % Initial temperature
    end

    function [pl, ql, pr, qr] = bcfun(xl, ul, xr, ur, t)
        pl = h_hat * (u_in - ul);
        ql = 1;
        pr = 0;
        qr = 1;
    end
end
