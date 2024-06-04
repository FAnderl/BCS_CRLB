
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
%   2) UPDATE: Numerical Calculations
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
syms transcription_translation(t,TF)                                % Transcription/ Activation       

assumeAlso([gamma(t) degrad_S_ext(t,S_ext) degrad_S_int(t,S_int) degrad_RR(t,RR) degrad_TF(t,TF) degrad_P(t,P) active_up(t,S_ext) basal_RR(t,RR) TF_maturation(t,S_int, RR) TF_dematuration(t,TF) transcription_translation(t,TF)], 'real')

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
gamma(t) = diff(gamma_const * exp(-((t-200)^2)/50),t);

% Placeholder Argument with the right dimensions
syms u [NUM_STATE_VARS 1] real;


% Vector a of Propensity Functions
a_vec(u) = [gamma(t); degrad_S_ext(t, u(1)); active_up(t, u(1)); degrad_S_int(t, u(2)); basal_RR; ...
    TF_maturation(t, u(2), u(3)); degrad_RR(t, u(3)); TF_dematuration(t, u(4));  degrad_TF(t, u(4)); transcription_translation(t, u(4)); degrad_P(t, u(5))];

a_vec_x = a_vec(x_vec(1),x_vec(2),x_vec(3), x_vec(4), x_vec(5));


% Covariance Matrix
syms Sigma real
Sigma = diag(S * a_vec(x_vec(1),x_vec(2),x_vec(3), x_vec(4), x_vec(5)));    % TODO: Unused -> Remove


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

delta_S = 0.003;
delta_RR = 0.01;
delta_TF = 0.1;
delta_P = 0.1;

kappa_up = 0.3;
K_up = 0.1;
kappa_RR = 0.1;
k_TF_m = 0.1;
k_TF_f = 0.1;
kappa_T = 0.5;
K_T = 0.1;
alpha_T = 1;



gamma_const = 100;


% Observation Noise Variance
sigma_obs = 1;
xi = 1;

% Define time step tau
tau = 1;


init_states = [1e-3, 1e-3, 100, 1e-3, 1e-3]

% Solve ODE system 

F_ode=@(t,Y) F(t, Y, eval(subs(gma_params)));
t0 = 0;
t_end = 500;
y0 = init_states;  % TODO: Replace with FLAG coded intital conditions; Is that a problem?
opt = odeset('RelTol',1e-6,'AbsTol',1e-6,'Vectorized','on',NonNegative=[1,2,3,4,5]);
figure(10)
gma_sol = ode15s(F_ode,[t0:tau:t_end    ],y0,opt);
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
S_ext = init_states(1)
S_int = init_states(2)
RR = init_states(3)
TF = init_states(4)
P = init_states(5)

% Initia time for J0 - > J1
t = 0

% J1
subs(D_11)
subs(D_21)
subs(D_12)
subs(D_22)



J = eval(D_22) - (eval(D_21)*pinv(J_0 + eval(D_11))) * eval(D_12);


% Container for FIM for each time step

J_FIM = zeros(NUM_STATE_VARS, NUM_STATE_VARS, 250);

inv_J_FIM =zeros(NUM_STATE_VARS, NUM_STATE_VARS, 250);

num_t_steps = size(sol_t,2);

for i = 2:num_t_steps

    % Update States - J_{n+1}
    S_ext = sol_t(1, i);
    S_int = sol_t(2,i);
    RR = sol_t(3,i);
    TF = sol_t(4,i);
    P = sol_t(5,i); 


    % Guard against *DIVISION BY ZERO*
    if S_ext == 0
        S_ext = 0.0001;
    end
    if S_int == 0
        S_int = 0.0001;
    end

    if RR == 0
        RR = 0.0001;  
    end
    if TF == 0
        TF = 0.0001;
    end
    if P == 0
        P = 0.0001;
    end
    


    % Set t (only relevant if gamma(t) is a function of t)
    time_vec = t0:tau:t_end;
    t = time_vec(i);



    % Substitute - J_{n+1}

    subs(D_11);
    subs(D_21);
    subs(D_12);
    subs(D_22);


    J = eval(D_22) - (eval(D_21)*pinv(J + eval(D_11))) * eval(D_12);

    J_FIM(:,:,i) = J;
    
    inv_J_FIM(:,:,i) = pinv(J);


end

%% Plotting

% Plot the element of the inverse FIM corresponding to S_ext (first element)
figure(11)
el = inv_J_FIM(1,1,:);
plot(1:num_t_steps, squeeze(el))
title("Inverse FIM Element (1,1) over Time")
xlabel("Time")
ylabel("Inverse FIM Element (1,1)")





%% Easy facilitated toy example with fewer states and simple propensity function vector

clear
% Define State Variables 
syms X1 X2 real
% Vector x of State Variables
x_vec = [X1; X2];
assumeAlso(x_vec,'real')

% Parameters 
syms k1 k2 d1 d2 real
params = [k1, k2, d1, d2];

% Stoichiometry Matrix
syms S real
S = sym([
    1 -1 -1 0 ;
    0 0  1 -1 
]);

% Define Propensity Functions 
syms a1(t, X1, X2) a2(t, X1, X2) real
a1(t, X1) = k1;
deg1(t, X1) = d1 * X1; 
a2(t, X1) = k2 * X1;
deg2(t, X2) = d2 * X2; 

% Placeholder Argument with the right dimensions
syms u [2 1] real;

% Vector a of Propensity Functions
a_vec(u) = [a1(t, u(1));deg1(t, u(1)) ; a2(t, u(1)); deg2(t, u(2))];


% Define State Space System
F(u) = u + S * a_vec(u(1), u(2));

% Mean 
z_np1 = F(x_vec(1), x_vec(2));

% Covariance Matrix
G(u) = diag(S * sqrt(a_vec(u(1), u(2))));
v_np1 = G(x_vec(1), x_vec(2));
Sigma_n = v_np1 * v_np1.';

% D^11
D_11 = sym(zeros(2, 2), 'r');
for i = 1:2
    for j = 1:2
        D_11(i,j) = diff(z_np1.', x_vec(i)) * pinv(Sigma_n) * diff(z_np1, x_vec(j))  + ...
            1/2 * trace(pinv(Sigma_n) * diff(Sigma_n,x_vec(i)) * pinv(Sigma_n)* diff(Sigma_n,x_vec(j)));
    end
end

% D^12
syms D_12  real % Matrix (2x2)
D_12 = - jacobian(z_np1.',x_vec) * pinv(Sigma_n) ;

% D^21
D_21 = D_12.';

% Observation Model (assumed Gaussian)
syms xi sigma_obs real

% D^22
h = [0; xi*X2];    % Observation vector
D_22 = inv(Sigma_n) +  (jacobian(h.', x_vec)  *   1/(sigma_obs^2)  * jacobian(h.', x_vec).') ;

% Setting up ODE state space model for solving
syms X1(t) X2(t)   % state variables
rhs = [diff(X1(t), t); diff(X2(t), t)] == S*a_vec(x_vec(1), x_vec(2));

% Defining Symbolic ODE System
diff_X1 = rhs(1);
diff_X2 = rhs(2);
eq_sys = [diff_X1; diff_X2];
sys_vars =["X1", "X2"];
vars = [X1(t), X2(t)];

[MS,F] = massMatrixForm(eq_sys,vars);
MS = odeFunction(MS,vars);
F = odeFunction(F,vars,params);

% Substitution with real, numerical values
% Substitute Parameter Values
k1 = 0.1;
k2 = 0.2;
d1 = 0.1; 
d2 = 0.1; 

% Observation Noise Variance
sigma_obs = 1;
xi = 1;

% Define time step tau
tau = 1;

init_states_toy = [1, 1e-4];

% Solve ODE system 
F_ode=@(t,Y) F(t, Y, eval(subs(params)));
t0 = 0;
t_end = 10;
y0 = init_states_toy;

opt = odeset('RelTol',1e-6,'AbsTol',1e-6,'Vectorized','on',NonNegative=[1,2]);
figure(10)
sol = ode15s(F_ode,[t0:tau:t_end],y0,opt);
ode15s(F_ode,[t0:tau:t_end],y0,opt);
legend('X1', 'X2')
title("Toy Example Dynamic Solution")


sol_toy = deval(t0:tau:t_end, sol);

% FIM RECURSION

% Define the (initial ?) states

J_0 = (zeros(2, 2));

% Initial values for the states for computation of J_1; Values cannot be
% zero due to *DIVISION BY ZERO*

X1 = init_states_toy(1);
X2 = init_states_toy(2);


% Initia time for J0 - > J1
t = 0;

% J1
subs(D_11);
subs(D_21);
subs(D_12);
subs(D_22);

J = eval(D_22) - (eval(D_21)*pinv(J_0 + eval(D_11))) * eval(D_12);

% Container for FIM for each time step
J_FIM = zeros(2, 2, 10);
inv_J_FIM =zeros(2, 2, 10);

num_t_steps = size(sol_toy,2);

for i = 2:num_t_steps

    % Update States - J_{n+1}
    X1 = sol_toy(1, i);
    X2 = sol_toy(2,i);

    % Guard against *DIVISION BY ZERO*
    if X1 == 0
        X1 = 0.0001;
    end
    if X2 == 0
        X2 = 0.0001;
    end

    % Set t (only relevant if gamma(t) is a function of t)
    time_vec = t0:tau:t_end;
    t = time_vec(i);

    % Substitute - J_{n+1}
    subs(D_11);
    subs(D_21);
    subs(D_12);
    subs(D_22);

    J = eval(D_22) - (eval(D_21)*pinv(J + eval(D_11))) * eval(D_12);

    J_FIM(:,:,i) = J;
    inv_J_FIM(:,:,i) = pinv(J);

end



% Plot the element of the inverse FIM corresponding to X1 (first element)
figure(11)
el = inv_J_FIM(2,2,:);
plot(1:num_t_steps, squeeze(el))
title("Inverse FIM Element (1,1) over Time")
xlabel("Time")
ylabel("Inverse FIM Element (1,1)")

