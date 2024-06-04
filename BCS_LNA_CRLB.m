
%%%%%%%%%%%%%%%%%%%%%%%%%%%% BCS_LNA_CRLB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Info: This script serves as code base for the planned IEEE Sensors 2024
%   conference paper
%
%   The main idea is to computer the recursive CRLB for the estimation of
%   a latent state in a state space model based on stochastic biochemical
%   kinetics
%
%   Specifically the following things are implemented here: 
%   1) Symbolic calculation of derivatives of log-likelihood function
%   2) UPDATE: Symbolic Calculationb 
%
%    Authors: Florian Anderl (florian.anderl@ntnu.no)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% GLOBAL FLAGS

NUM_STATE_VARS = 5              % UNUSED at the mooment :(
NUM_REACTION_CHANNELS = 11       % loosely referring to the protein transduction model

% time 
syms t tau real


 % Define State Variables 
syms S_ext S_int RR TF P real


% Vector x of State Variables
x_vec = [S_ext; S_int; RR; TF; P];

assumeAlso(x_vec,'real')

% Parameters 
syms delta_S   real
syms delta_RR  real
syms delta_TF  real
syms delta_P   real

syms kappa_up  real
syms K_up      real
syms kappa_RR  real
syms k_TF_m    real
syms k_TF_f    real
syms kappa_T   real
syms K_T       real
syms alpha_T   real       % Translation Efficiency

syms gamma_const real


gma_params = [delta_S, delta_RR, delta_TF, delta_P, kappa_up, K_up, kappa_RR, k_TF_m, k_TF_f, kappa_T, K_T, alpha_T, gamma_const];


% Stoichiometry Matrix
syms S real
S = sym([
    1 -1 -1 0 0 0 0 0 0 0 0;
    0 0 1 -1 0 -2 0 2 0 0 0;
    0 0 0 0 1 -2 -1 2 0 0 0;
    0 0 0 0 0 1 0 -1 -1 0 0;
    0 0 0 0 0 0 0 0 0 1 -1
]);


% Define Propensity Functions 
syms gamma(t)                                                       % Influx of S_ext
syms degrad_S_ext(t,S_ext)                                          % Degradation of S_ext
syms degrad_S_int(t,S_int)                                          % Degradation of S_int
syms degrad_RR(t,RR)                                                % Degradation of RR
syms degrad_TF(t,TF)                                                % Degradation of TF
syms degrad_P(t,P)                                                  % Degradation of P

syms active_up(t,S_ext)                                             % Active Uptake of S_ext (Michaelis-Menten)
syms basal_RR(t,RR)                                                 % Basal Production of RR
syms TF_maturation(t,S_int, RR)                                     % Maturation of TF
syms TF_dematuration(t,TF)                                          % Degradation of TF
syms transcription_translation(t,TF)                                                % Transcription/ Activation       

assumeAlso([gamma(t) degrad_S_ext(t,S_ext) degrad_S_int(t,S_int) degrad_RR(t,RR) degrad_TF(t,TF) degrad_P(t,P) active_up(t,S_ext) basal_RR(t,RR) TF_maturation(t,S_int, RR) TF_dematuration(t,TF) active_TF(t,TF)], 'real')

degrad_S_ext(t,S_ext) = delta_S * S_ext;
degrad_S_int(t,S_int) = delta_S * S_int;
degrad_RR(t,RR) = delta_RR * RR;
degrad_TF(t,TF) = delta_TF * TF;
degrad_P(t,P) = delta_P * P;
active_up(t,S_ext) = kappa_up * S_ext / (K_up + S_ext);
basal_RR(t,RR) = kappa_RR;
TF_maturation(t,S_int, RR) = k_TF_m * S_int * RR;
TF_dematuration(t,TF) = k_TF_f * TF;
transcription_translation(t,TF) = alpha_T * kappa_T * TF / (K_T + TF);
gamma(t) = gamma_const;

% Placeholder Argument with the right dimensions
syms u [NUM_STATE_VARS 1] real;


% Vector a of Propensity Functions
a_vec(u) = [gamma(t); degrad_S_ext(t, u(1)); active_up(t, u(1)); degrad_S_int(t, u(2)); basal_RR; ...
    TF_maturation(t, u(2), u(3)); degrad_RR(t, u(3)); TF_dematuration(t, u(4));  degrad_TF(t, u(4)); transcription_translation(t, u(4)); degrad_P(t, u(5))];

a_vec_x = a_vec(x_vec(1),x_vec(2),x_vec(3), x_vec(4), x_vec(5));


% Covariance Matrix
syms Sigma real
Sigma = diag(S * a_vec(x_vec(1),x_vec(2),x_vec(3), x_vec(4), x_vec(5)));


% Define State Space System

% TODO: I need to establish 
F(u) = u + S * a_vec(u(1), u(2), u(3), u(4), u(5)) * tau;

% Mean 
z_np1 = F(x_vec(1),x_vec(2),x_vec(3), x_vec(4), x_vec(5));

% Covariance Matrix
G(u) = diag(S * sqrt(a_vec(u(1), u(2), u(3), u(4), u(5)))) * tau;

v_np1 = G(x_vec(1),x_vec(2),x_vec(3), x_vec(4), x_vec(5));

Sigma_n = v_np1 * v_np1.';




% D^11

D_11 = sym(zeros(NUM_STATE_VARS, NUM_STATE_VARS), 'r');
for i = 1:NUM_STATE_VARS
    for j = 1:NUM_STATE_VARS
    D_11(i,j) = diff(z_np1.', x_vec(i)) * pinv(Sigma_n) * diff(z_np1, x_vec(j))  + ...
        1/2 * trace(pinv(Sigma_n) * diff(Sigma_n,x_vec(i)) * pinv(Sigma_n)* diff(Sigma_n,x_vec(j)));

    end
end


% D^12

syms D_12  real % Matrix (nxn)

D_12 = - jacobian(z_np1.',x_vec) * pinv(Sigma_n) ;





% D^21
D_21 = D_12.';





% Observation Model (assumed Gaussian)
syms xi sigma_obs real
 

% D^22
h = [0;0;0;0;xi*P]    % Observation vector



D_22 = inv(Sigma_n) +  (jacobian(h.', x_vec)  *   1/(sigma_obs^2)  * jacobian(h.', x_vec).') ;






% Setting up ODE state space model for solving

syms S_ext(t)  S_int(t) RR(t) TF(t) P(t)   % state variables
rhs = [diff(S_ext(t), t); diff(S_int(t), t); diff(RR(t), t); diff(TF(t), t); diff(P(t), t)] == S*a_vec(x_vec(1),x_vec(2),x_vec(3), x_vec(4), x_vec(5));


% Defining Symbolic ODE System
diff_S_ext = rhs(1)
diff_S_int = rhs(2)
diff_RR = rhs(3)
diff_TF = rhs(4)
diff_P  = rhs(5)

eq_sys = [diff_S_ext; diff_S_int; diff_RR; diff_TF; diff_P];


sys_vars =["S_ext","S_int", "RR", "TF", "P"];

vars = [S_ext(t),S_int(t), RR(t), TF(t), P(t)]; % TODO: Fix naming wrt. sys_vars and vars




[MS,F] = massMatrixForm(eq_sys,vars);

MS = odeFunction(MS,vars);

F = odeFunction(F,vars,gma_params);





% Substitution with real, numerical values

% Substitute Parameter Values

delta_S = 0.1;
delta_RR = 0.1;
delta_TF = 0.1;
delta_P = 0.1;

kappa_up = 0.1;
K_up = 0.1;
kappa_RR = 0.1;
k_TF_m = 0.1;
k_TF_f = 0.1;
kappa_T = 0.5;
K_T = 0.1;
alpha_T = 1;

gamma_const = 0.1;


% Observation Noise Variance
sigma_obs = 1;
xi = 1;

% Define time step tau
tau = 0.1;


% Solve ODE system 

F_ode=@(t,Y) F(t, Y, eval(subs(gma_params)));
t0 = 0;
t_end = 500;
y0 = [100,1,100,1,1];  % TODO: Replace with FLAG coded intital conditions; Is that a problem?
opt = odeset('RelTol',1e-6,'AbsTol',1e-6,'Vectorized','on');
figure(10)
gma_sol = ode15s(F_ode,[t0:tau:t_end],y0,opt);
ode15s(F_ode,[t0:tau:t_end],y0,opt);
legend('S_ext','S_int', 'RR', 'TF', 'P')
title("Original GMA-system Dynamic Solution")



% Evaluate solution at the desired intervals (tau)
sol_t = deval(t0:tau:t_end, gma_sol)





%% FIM RECURSION

% Define the (initial ?) states
% TODO: How to choose  J_0?
J_0 = (zeros(NUM_STATE_VARS, NUM_STATE_VARS));

% Initial values for the states for computation of J_1; Values cannot be
% zero due to *DIVISION BY ZERO*
S_ext = 100
S_int = 1
RR = 100
TF = 1
P = 1

% J1
subs(D_11)
subs(D_21)
subs(D_12)
subs(D_22)

J = eval(D_22) - (eval(D_21)*pinv(J_0 + eval(D_11))) * eval(D_12);


% Container for FIM for each time step

J_FIM = zeros(NUM_STATE_VARS, NUM_STATE_VARS, 250);

inv_J_FIM =zeros(NUM_STATE_VARS, NUM_STATE_VARS, 250);

for i = 2:250

    % Update States - J_{n+1}
    S_ext = sol_t(1, i);
    S_int = sol_t(2,i);
    RR = sol_t(3,i);
    TF = sol_t(4,i);
    P = sol_t(5,i); 

    % Substitute - J_{n+1}

    subs(D_11);
    subs(D_21);
    subs(D_12);
    subs(D_22);


    J = eval(D_22) - (eval(D_21)*pinv(J + eval(D_11))) * eval(D_12);

    J_FIM(:,:,i) = J;
    
    inv_J_FIM(:,:,i) = pinv(J);


end

%% Ploting

% Plot the element of the inverse FIM corresponding to S_ext (first element)
figure(11)
el = inv_J_FIM(3,3,:);
plot(1:250, squeeze(el))
title("Inverse FIM Element (1,1) over Time")
xlabel("Time")
ylabel("Inverse FIM Element (1,1)")




