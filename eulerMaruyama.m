function [t, X] = eulerMaruyama(F_drift, G_diffusion, init_states, num_stoch ,dt, t_end)
    % eulerMaruyama Simulates an SDE using the Euler-Maruyama method
    % 
    % Inputs:
    %   F_drift - Function handle for the drift term (F)
    %   G_diffusion - Function handle for the diffusion term (G)
    %   init_states - Initial states of the system
    %   dt - Time step size
    %   t_end - End time of the simulation
    %
    % Outputs:
    %   t - Time vector
    %   X - Matrix of state values over time

    % Number of steps
    numSteps = floor(t_end / dt);
    % Number of state variables
    numVars = length(init_states);
    
    % Preallocate arrays for efficiency
    t = (0:dt:t_end)';
    X = zeros(numSteps + 1, numVars);
    
    
    % Set initial state
    X(1, :) = init_states;
    
    % Simulation loop
    for i = 1:numSteps
        % Get current state
        Xt = X(i, :)';
        % Calculate drift and diffusion
        drift = F_drift(t(i), Xt);
        diffusion = G_diffusion(t(i), Xt);
        % Generate Wiener process increments
        dW = sqrt(dt) * randn(num_stoch, 1);
        % Update state with Euler-Maruyama method
        X(i+1, :) = Xt + drift * dt + diffusion * dW;
        % Enforce non-negativity
        X(i+1, :) = max(X(i+1, :), 0);
    end
end
