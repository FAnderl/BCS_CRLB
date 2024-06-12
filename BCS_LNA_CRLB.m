
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
%   3) Ultimately, this         
%
%    Authors: Florian Anderl (florian.anderl@ntnu.no)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main
clear
%GLOBAL FLAGS

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
% TODO: 1,2 : 0 -> -1 
S = sym([
    1 -1 -1 0 0 0 0 0 0 0 0;
    0 0 1 -1 0 -2 0 2 0 0 0;
    0 0 0 0 1 -2 -1 2 0 0 0;
    0 0 0 0 0 1 0 -1 -1 0 0;
    0 0 0 0 0 0 0 0 0 1 -1
])


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

degrad_S_ext(t,S_ext) =  delta_S * S_ext;
degrad_S_int(t,S_int) = delta_S * S_int;
degrad_RR(t,RR) = delta_RR * RR;
degrad_TF(t,TF) = delta_TF * TF;
degrad_P(t,P) = delta_P * P;
active_up(t,S_ext) = kappa_up * S_ext / (K_up + S_ext);
basal_RR(t,RR) = kappa_RR;
TF_maturation(t,S_int, RR) = k_TF_m * S_int^2 * RR^2;
TF_dematuration(t,TF) = k_TF_f * TF;
transcription_translation(t,TF) = alpha_T * kappa_T * TF / (K_T + TF);
gamma(t) = diff(gamma_const/2 * sin(2*pi*0.1*t));%diff(gamma_const * exp(-((t-200)^2)/200),t);%; %diff(gamma_const * exp(-((t-200)^2)/50),t);

% Placeholder Argument with the right dimensions
syms u [NUM_STATE_VARS 1] real;


% Vector a of Propensity Functions
a_vec(u) = [gamma(t); degrad_S_ext(t, u(1)); active_up(t, u(1)); degrad_S_int(t, u(2)); basal_RR; ...
    TF_maturation(t, u(2), u(3)); degrad_RR(t, u(3)); TF_dematuration(t, u(4));  degrad_TF(t, u(4)); transcription_translation(t, u(4)); degrad_P(t, u(5))];

a_vec_x = a_vec(x_vec(1),x_vec(2),x_vec(3), x_vec(4), x_vec(5));





% Covariance Matrix
syms Sigma real


% Define State Space System

% Mean State Equation
F(u) = u + S * a_vec(u(1), u(2), u(3), u(4), u(5)) * tau;

% Mean 
z_np1 = F(x_vec(1),x_vec(2),x_vec(3), x_vec(4), x_vec(5));

% Covariance Matrix
G(u) = diag(S * sqrt(a_vec(u(1), u(2), u(3), u(4), u(5)))) * tau;


% This is probably obsoltet now; remove consequently...
v_np1 = G(x_vec(1),x_vec(2),x_vec(3), x_vec(4), x_vec(5));
Sigma_n = v_np1 * v_np1.';



% Redefine Sigma according to updated model
Sigma = sym(zeros(NUM_STATE_VARS, NUM_STATE_VARS),'r');

for i = 1 : NUM_STATE_VARS
    for j = 1 : NUM_STATE_VARS
   
        % sigma_ii
            if i == j
                  temp = S.^2*a_vec_x; 
                  % Could easily be described by S*a (this should also be done in this fashion in the paper...done!)
                  Sigma(i,j) =  tau * temp(i);
                
                
            elseif i ~= j 
                    

                    for k = 1 : NUM_REACTION_CHANNELS

                        Sigma(i,j) = Sigma(i,j) + (tau * (S(i,k)* S(j,k)) * (a_vec_x(k)));
                  
                    end
            
            end
    end
end


% This is due to implementation-specific reasons; ultimately make this
% nicer...
Sigma_n = Sigma

% If all fails, make Sigma Idendity Matrix
% Sigma_n = sym(eye(NUM_STATE_VARS,NUM_STATE_VARS));

% Inverse of Sigma % THIS MIGHT BE PROBLEMATIC
i_Sigma_n = inv(Sigma_n);

% Define a place holder matrix i_Sigma_n and then calculate the inversion
% ONLY NUMERICALLY
syms i_Sigma_n [NUM_STATE_VARS,NUM_STATE_VARS]

% D^11

D_11 = sym(zeros(NUM_STATE_VARS, NUM_STATE_VARS), 'r');
for i = 1:NUM_STATE_VARS
    for j = 1:NUM_STATE_VARS
    D_11(i,j) = diff(z_np1, x_vec(i)).' * i_Sigma_n * diff(z_np1, x_vec(j))  + ...
        (1/2) * trace(i_Sigma_n * diff(Sigma_n,x_vec(i)) * i_Sigma_n* diff(Sigma_n,x_vec(j)));

    end
end


% D^12

syms D_12  real % Matrix (nxn)

% THIS has been changed. I think I have to transpose the jacobian directly!
D_12 = - jacobian(z_np1,x_vec).' * i_Sigma_n ;





% D^21
D_21 = D_12.';





% Observation Model (assumed Gaussian)
syms xi sigma_obs real
 

% D^22
h = [0;0;0;0;xi*P]    % Observation vector



D_22 = i_Sigma_n +  (jacobian(h.', x_vec).'  *   1/(sigma_obs^2)  * jacobian(h.', x_vec)) ;






%%  Setting up ODE state space model for solving & Solving it

syms S_ext(t)  S_int(t) RR(t) TF(t) P(t)   % state variables
rhs = [diff(S_ext(t), t); diff(S_int(t), t); diff(RR(t), t); diff(TF(t), t); diff(P(t), t)] == S*a_vec(x_vec(1),x_vec(2),x_vec(3), x_vec(4), x_vec(5));


% Defining Symbolic ODE System
diff_S_ext = rhs(1)
diff_S_int = rhs(2)
diff_RR = rhs(3)
diff_TF = rhs(4)
diff_P  = rhs(5)

eq_sys = [diff_S_ext; diff_S_int; diff_RR; diff_TF; diff_P];


sys_vars =["S_{ext}","S_{int}", "RR", "TF", "P"];

vars = [S_ext(t),S_int(t), RR(t), TF(t), P(t)]; % TODO: Fix naming wrt. sys_vars and vars




[MS,F] = massMatrixForm(eq_sys,vars);

MS = odeFunction(MS,vars);

F = odeFunction(F,vars,gma_params);





% Substitution with real, numerical values

% Substitute Parameter Values

delta_S = 0.003;
delta_RR = 0.01;
delta_TF = 0.1;
delta_P = 0.02;

kappa_up = 1;
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


init_states = [gamma_const, 1, 100, 1, 1];

% Solve ODE system 

F_ode=@(t,Y) F(t, Y, eval(subs(gma_params)));
t0 = 0;
t_end = 600;
y0 = init_states;  % TODO: Replace with FLAG coded intital conditions; Is that a problem?
opt = odeset('RelTol',1e-6,'AbsTol',1e-6,'Vectorized','on',NonNegative=[1,2,3,4,5],OutputFcn='odeplot');
figure(10)
gma_sol = ode15s(F_ode,[t0:tau:t_end],y0,opt);
%ode15s(F_ode,[t0:tau:t_end],y0,opt);
legend('S_ext','S_int', 'RR', 'TF', 'P')
title("Original GMA-system Dynamic Solution")



% Evaluate solution at the desired intervals (tau)
sol_t = deval(t0:tau:t_end, gma_sol);





%% FIM RECURSION

% Define the (initial ?) states
% TODO: How to choose  J_0?
J_0 = (zeros(NUM_STATE_VARS, NUM_STATE_VARS));

% Initial values for the states for computation of J_1; Values cannot be
% zero due to *DIVISION BY ZERO*
S_ext = init_states(1);
S_int = init_states(2);
RR = init_states(3);
TF = init_states(4);
P = init_states(5);

% Initial time for J0 - > J1
t = 0;



num_i_Sigma_n = pinv(eval(subs(Sigma_n)));

% Substitution of i_Sigma_n 
for i = 1: NUM_STATE_VARS
    for j = 1: NUM_STATE_VARS
        eval(strcat(string(i_Sigma_n(i,j)),"=",string(num_i_Sigma_n(i,j))));
    end
end


% J1
subs(D_11)
subs(D_21)
subs(D_12)
subs(D_22)


% Tichavsky Recursion Equation
J = eval(D_22) - (eval(D_21)*pinv(J_0 + eval(D_11))) * eval(D_12);


% Container for FIM for each time step

J_FIM = zeros(NUM_STATE_VARS, NUM_STATE_VARS, 250);

inv_J_FIM =zeros(NUM_STATE_VARS, NUM_STATE_VARS, 250);

num_t_steps = size(sol_t,2);


% Some restructuring so that I can use an optimized funtion for recursion
% computation
D_11_cont = zeros(NUM_STATE_VARS, NUM_STATE_VARS, num_t_steps);
D_12_cont = zeros(NUM_STATE_VARS, NUM_STATE_VARS, num_t_steps);
D_21_cont = zeros(NUM_STATE_VARS, NUM_STATE_VARS, num_t_steps);
D_22_cont = zeros(NUM_STATE_VARS, NUM_STATE_VARS, num_t_steps);

% Fill up the D-containers
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


    num_i_Sigma_n = pinv(eval(subs(Sigma_n)));

    % Substitution of i_Sigma_n 
    for k = 1: NUM_STATE_VARS
        for j = 1: NUM_STATE_VARS
            eval(strcat(string(i_Sigma_n(k,j)),"=",string(num_i_Sigma_n(k,j)),";"));
        end
    end


    D_11_cont(:,:,i) = eval(subs(D_11));
    D_21_cont(:,:,i) = eval(subs(D_21));
    D_12_cont(:,:,i) = eval(subs(D_12));
    D_22_cont(:,:,i) = eval(subs(D_22));


end

% This creates the optimized MEX function (It is a bit suboptimal that this is done everytime atm but accounts for the ARGS having possibly different dimensions depending on tau; should be optimized in the future)
% eval("codegen calculateRecursion.m -args {D_11_cont,D_12_cont,D_21_cont,D_22_cont, J,NUM_STATE_VARS} ")

%Calculate Recursion
[J_FIM, J_FIM_inv]  = calculateRecursion_mex(D_11_cont, D_12_cont, D_21_cont, D_22_cont, J, NUM_STATE_VARS); 

% Calculate the recursion 
%[J_FIM, J_FIM_inv] = calculateRecursion(D_11_cont, D_12_cont, D_21_cont, D_22_cont, J, NUM_STATE_VARS);


%% Plotting

% Plot the element of the inverse FIM corresponding to S_ext (first element)
figure(11)
for i = 1:NUM_STATE_VARS
    subplot(2,3,i)
    el = J_FIM_inv(i,i,:);
    plot(1:num_t_steps, squeeze(el), 'LineWidth', 2)
    title(sprintf("Inverse FIM Element (%d,%d) over Time", i, i))
    xlabel("Time")
    ylabel(sprintf("Inverse FIM Element (%d,%d)", i, i))
    grid on
    set(gca, 'LineWidth', 1.5)
end
sgtitle("Inverse FIM Elements over Time")

% Log-pots
figure(12)
for i = 1:NUM_STATE_VARS
    subplot(2,3,i)
    el = J_FIM_inv(i,i,:);
    plot(1:num_t_steps, log(squeeze(el)), 'LineWidth', 2)
    title(sprintf("Log-Inverse FIM Element (%d,%d) over Time", i, i))
    xlabel("Time")
    ylabel(sprintf("Log-Inverse FIM Element (%d,%d)", i, i))
    grid on
    set(gca, 'LineWidth', 1.5)
end