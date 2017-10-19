function d = usadel2(u)
    % Usage: give the program the phase-winding u=ξ(∂φ/∂z) as
    % its input, and it will then self-consistently solve the
    % Usadel equation in a bulk superconductor with this u.
    % It uses this calculation to plot the density of states.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Predefine constants and arrays
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Energies ϵ/Δ₀
    E = [linspace(1e-3,  1.5, 400), ...
         linspace(1.501, 4.5, 100), ...
         linspace(4.501, 30,  100)] + 0.01i;

    % Solutions Φ(ϵ)
    F = zeros(size(E));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform the calculation itself
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initial guess for the gap δ
    d  = 1;

    % Storage for previous results
    d_ = 0;

    % Self-consistency routine
    while abs(d - d_) > 1e-4
        % Solve the fixpoint iteration
        for n=length(E):-1:1
            % Initial guess for Φ(ϵ)
            if n == length(E)
                F(n) = 1;
            else
                F(n) = F(n+1);
            end

            % Storage for previous results
            F_ = 0;

            % Newton's method for finding the root
            while abs(F(n) - F_) > 1e-5
                F_   = F(n);
                F(n) = F(n) - f(F(n), E(n), d, u)/df(F(n), E(n), d, u);
            end
        end

        % Calculate the gap integrand
        i = F./sqrt(E.^2 - F.^2);

        % Perform the selfconsistency integral the "naive way"
        I = sum( real(E(2:end)-E(1:end-1)) .* real(i(2:end)+i(1:end-1)) )/2;

        % Use this to update the gap estimate
        d_ = d;
        d  = I/acosh(real(E(end)));

        % Status information
        disp(['Current gap: ', num2str(d)]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Post-processing of results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calculate the density of states
    D = real(E./sqrt(E.^2 - F.^2));
    
    % Visualize the final results
    plot([real(fliplr(-E)) real(E)], [fliplr(D) D]);
    xlim([-2, 2]);
    ylim([ 0, 3]);
    xlabel('Energy \epsilon/\Delta_0');
    ylabel('Density of states D(\epsilon)');
end

function y = f(F, E, d, u)
    % This function defines the fixpoint function
    % f(Φ) = 0 that we are looking for roots of.
    y = F - d/(1+u^2/(2*sqrt(F^2 - E^2)));
end

function y = df(F, E, d, u)
    % This function makes a crude estimate of the
    % derivative of the fixpoint function above.
    dF = 1e-6;;
    y = (f(F+dF/2, E, d, u) - f(F-dF/2, E, d, u))/dF;
end
