# Packet-Loss-Delay-Analysis-Wireless-Fading-Channel
packet loss and end-to-end delay analysis for wireless fading channels using Stochastic Network Calculus (SNC) framework

* **Project Scope:** High-level summary of how the paper analyzes Quality of Service (QoS) constraints for Ultra-Reliable Low-Latency Communication (URLLC) networks.
* **Core Methodology:** Clean mathematical formulas detailing the mapping of a Rayleigh fading channel into a Bernoulli ON-OFF channel, along with its throughput optimization function.
* **Network Calculus Bounds:** Formulations for the Stochastic Arrival Curve (Compound Poisson process) and Weak Stochastic Service Curve, alongside the exact equations derived for **Packet Loss** and **Delay Distribution**.
* **Key Insights:** Brief summary of performance sensitivities regarding the delayed-start time ($T$) and free tuning parameters ($z_\alpha, z_\beta$).

## code

% Packet Loss and Delay Analysis for a Wireless Fading Channel
clear; clc; close all;

% Global Parameters
W_base = 15e3;     
delta2 = 1;        
p      = 0.52;     
lambda_base = 2;   
l      = 1e3;      

% Compute R_opt
R_vec = linspace(0, 8e4, 5000); 
throughput = R_vec .* exp( -(2.^(R_vec ./ W_base) - 1) ./ (2 * delta2) );
[max_throughput, idx_max] = max(throughput);
R_opt_val = R_vec(idx_max); 

% Figure 1: Throughput vs Rate
figure('Name', 'Fig 1: Average Throughput', 'Color', 'w');
plot(R_vec, throughput, 'r', 'LineWidth', 2); hold on;
yline(max_throughput, 'k--', 'LineWidth', 1);
xlabel('Transmission rate (bps)');
ylabel('Average throughput (bps)');
title('Fig. 1. Relationship between average throughput and transmission rate');
grid on;
axis([0 8e4 0 10000]);

% Figure 2: Binary Function
figure('Name', 'Fig 2: Binary Function', 'Color', 'w');
za_vec = linspace(1e-5, 6e-4, 50);
zb_vec = linspace(1e-5, 6e-4, 50);
[ZA, ZB] = meshgrid(za_vec, zb_vec);

BinFunc_Val = (lambda_base ./ ZA) .* (exp(ZA .* l) - 1) + ...
              log(1 - p + p .* exp(-ZB .* R_opt_val)) ./ ZB;

surf(ZA, ZB, BinFunc_Val);
shading interp;
colormap jet;
xlabel('z_\alpha');
ylabel('z_\beta');
zlabel('Binary function');
title('Fig. 2. Numerical distribution of the binary function');
view(-45, 30);
grid on;

% Figure 3: Packet Loss Distribution
figure('Name', 'Fig 3: Packet Loss', 'Color', 'w');
hold on;

lambda_sim = 2000; 
za_loss = 2e-4; 
zb_loss = 2e-4;
T_list = [1e-3, 5e-3, 9e-3];
L_vec = linspace(0, 1e5, 200);

styles = {'-r', '-+b', '-*g'};

for i = 1:length(T_list)
    T_curr = T_list(i);
    burst_term = (lambda_sim / za_loss) * (exp(za_loss * l) - 1) * T_curr;
    arg_val = L_vec - burst_term;    
    P_loss = calculate_f_conv_g(arg_val, za_loss, zb_loss);
    plot(L_vec, P_loss, styles{i}, 'LineWidth', 1.5, 'MarkerIndices', 1:10:length(L_vec));
end

xlabel('Buffer size (bits)');
ylabel('Packet loss distribution');
title('Fig. 3. Packet loss distribution for different delayed-start times');
legend('1 ms','5 ms','9 ms');
grid on;
ylim([0 1]);

% Figure 4: Delay Distribution
figure('Name', 'Fig 4: Delay Distribution', 'Color', 'w');
hold on;

T_fig4 = 1;                  
t0_vec = linspace(0, 25, 60); 
z_cases = [1.5e-4, 2.0e-4, 3.0e-4];
fig4_styles = {'-r', '-+b', '-*g'};

for i = 1:length(z_cases)
    z_val = z_cases(i);
    psi = log(1 - p + p * exp(-z_val * R_opt_val));
    arg_delay = ((T_fig4 - t0_vec) ./ z_val) .* psi;
    P_delay = calculate_f_conv_g(arg_delay, z_val, z_val);
    plot(t0_vec, P_delay, fig4_styles{i}, 'LineWidth', 1.5, 'MarkerIndices', 1:3:length(t0_vec));
end

xlabel('Delay requirement (ms)');
ylabel('Delay distribution');
title('Fig. 4. Delay distribution for different z-values');
legend('1.5e-4','2.0e-4','3.0e-4');
grid on;
axis([0 25 0 1]);

% Figure 5: Delay Distribution Surface
figure('Name', 'Fig 5: Delay Surface', 'Color', 'w');

T_fig5 = 5;   
t0_fig5 = 10;
z_range = linspace(1e-5, 4e-4, 40);
[ZA_s, ZB_s] = meshgrid(z_range, z_range);

PSI_s = log(1 - p + p .* exp(-ZB_s .* R_opt_val));
ARG_s = ((T_fig5 - t0_fig5) ./ ZB_s) .* PSI_s;

denom_s = ZA_s + ZB_s;
log_term_s = log(ZA_s ./ ZB_s);

term1 = exp( -(ZA_s .* ZB_s .* ARG_s + ZA_s .* log_term_s) ./ denom_s );
term2 = exp( -(ZA_s .* ZB_s .* ARG_s - ZB_s .* log_term_s) ./ denom_s );

P_surf = term1 + term2;
P_surf(P_surf > 1) = 1;

surf(ZA_s, ZB_s, P_surf);
shading interp;
colormap jet;
xlabel('z_\alpha');
ylabel('z_\beta');
zlabel('Delay distribution');
title('Fig. 5. Delay distribution surface');
view(-45, 30);
grid on;

% Helper Function
function P = calculate_f_conv_g(x, za, zb)
    if za == zb
        log_val = 0;
    else
        log_val = log(za / zb);
    end
    
    denom = za + zb;

    E1 = -(za * zb .* x + za * log_val) ./ denom;
    E2 = -(za * zb .* x - zb * log_val) ./ denom;

    P = exp(E1) + exp(E2);
    P(P > 1) = 1;
end
