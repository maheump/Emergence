% This script finds finds the correspondance between the p(error) parameter
% (a probability substitution parameter used in Emergence) to the omega (w)
% (an exponential leak strength parameter used in previous projects). Both
% implements (exponentially) leaky integrations.
% See:
%   - Meyniel, F., Maheu, M., & Dehaene, S. (2016). Human inferences about
%       sequences: A minimal transition probability model. PLoS
%       computational biology, 12(12), e1005260.
%   - Maheu, M., Dehaene, S., & Meyniel, F. (2019). Brain signatures of a
%       multiscale process of sequence learning in humans. eLife, 8, e41541.
% 
% Copyright (c) 2020 Maxime Maheu

% Define the number of observations
N = 50;

% Define the grid for p(error)
ng = 1e3;
perrGrid = linspace(0,1/2,ng);

% Define the number of w parameters to test
maxomega = 50;

% Define the leak function
expdecayfun = @(w) (exp(-(1/w)*(N - (1:N)))) ./ 2 + 1/2;

% Prepare output variables
MSE = NaN(ng,maxomega);
bi  = NaN(1, maxomega);

% Over a grid of omega values
for omega = 1:maxomega
    
    % Get weights from the exponential decay function based on the current
    % value of omega
    leak = expdecayfun(omega);
    
    % Loop over the grid of p(error) parameter and compute corresponding
    % wieghts
    substitution = NaN(ng,N);
    for ig = 1:ng
        substitution(ig,:) = Emergence_IO_Leak(perrGrid(ig), N);
    end
    
    % Compute mean squared error between these two observation weights
    MSE(:,omega) = mean((substitution - leak) .^ 2, 2);
    
    % Find the value of p(error) that induce observation weights the
    % closest to the current value of omega
    [~,bi(omega)] = min(MSE(:,omega));
    
    % Display the results of the fit for one example w
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if omega == 16
        
        % Prepare a new window
        figure('Position', [1 955 250 150]);
        
        % Display the grid
        plot(1:N, substitution, '-', 'Color', g,  'LineWidth', 1/2); hold('on');
        
        % Display the lea
        plot(1:N, leak, 'k-', 'LineWidth', 2); hold('on');
        
        % Display the best fit to that leak
        plot(1:N, substitution(bi(omega),:), 'r--', 'LineWidth', 2); hold('on');
        
        % Customize the axes
        set(gca, 'TickLabelInterpreter', 'LaTeX');
        axis([1,N,1/2,1]);
        
        % Add some text labels
        xlabel('Observation $\#$', 'Interpreter', 'LaTeX');
        ylabel('$p(\mathrm{error})$', 'Interpreter', 'LaTeX');
        title(['$\omega = ', sprintf('%1.0f', omega), '\Leftrightarrow p(\mathrm{error}) = ', ...
            sprintf('%1.5f', perrGrid(bi(omega))), '$'], 'Interpreter', 'LaTeX');
        
    end
end

% Display w as a function of p(error)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare a new window
figure('Position', [252 955 250 150]);

% Display the relationship between leak omega parameter and the p(error)
% parameter
plot(1:maxomega, perrGrid(bi), 'k-', 'LineWidth', 2); hold('on');
scatter(1:maxomega, perrGrid(bi), 10, hsv(maxomega), 'filled');

% Customize the axes
set(gca, 'Box', 'Off', 'TickLabelInterpreter', 'LaTeX', 'XLim', [1,maxomega]);

% Add some text labels
xlabel('Leak parameter $\omega$', 'Interpreter', 'LaTeX');
ylabel('$p(\mathrm{error})$', 'Interpreter', 'LaTeX');
