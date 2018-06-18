% This script performs subject-specific analyses. In particular it produces
% a subject-specific figure specifying the trajectory of the finger in each
% sequence (both in the triangle and as a cumulative Barycentric
% coordinates plot) as well as the responses (s)he gave at the
% post-sequence questions.
% 
% Copyright (c) 2018 Maxime Maheu

% For each subject
for sub = 1:nSub
    fprintf('- Subject %2.0f/%2.0f\n', sub, nSub);
    F = cell(1,3);

    % For each type or regularity
    for reg = 1:3

        % Create an invisible figure
        f = figure('Visible', 'Off', 'Color', ones(1,3), 'Units', ...
            'Pixels', 'Position', [1 800 1920 400]);

        % For each sequence
        for seq = 1:numel(cidx{reg})
            cond = cidx{reg}(seq);

            % Display the triangle
            % ~~~~~~~~~~~~~~~~~~~~
            subplot(3, max(cellfun(@numel, cidx)), seq);

            % Display the trajectory on the triangle
            Emergence_PlotTrajOnTri(G{cond,sub}.BarycCoord, ... % trajectory
                                    floor(G{cond,sub}.Jump), ... % change point position
                                    tricol); % colors to use in the triangle

            % Display the barycentric coordinates
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            % Plot cumulative barycentric coordinates
            m = max(cellfun(@numel, cidx));
            subplot(3, m, seq + m);
            Emergence_PlotBarycTraj(G{cond,sub}.BarycCoord, tricol);

            % Display the change point position
            plot(repmat(G{cond,sub}.Jump,1,2), [0,1], 'k-', 'LineWidth', 1);

            % Customize the axes
            axis([1,N,0,1]);
            set(gca, 'YTick', 0:1/3:3, 'YTickLabel', {'0', '1/3', '2/3', '1'});
            set(gca, 'LineWidth', 2);

            % Display the inferred change point and detection point positions
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            % Display the inferred change point position
            subplot(3, m, seq + m*2);

            % Display the inferred change point position
            M = G{cond,sub}.Questions(3);
            E = (((100-G{cond,sub}.Questions(4)) / 100) * S{sub}.Nstims) / 2;
            plot(M, 1/2, 'ko', 'MarkerFaceColor', 'k'); hold('on'); hold('on');
            plot([M-E,M+E], repmat(1/2,1,2), 'k-', 'LineWidth', 1);

            % Display the inferred detection point position
            M = G{cond,sub}.Questions(5);
            E = ((100-G{cond,sub}.Questions(6)) / 100  *S{sub}.Nstims)/2;
            plot(M, 1/4, 'ko', 'MarkerFaceColor', 'k'); hold('on'); hold('on');
            plot([M-E,M+E], repmat(1/4,1,2), 'k-', 'LineWidth', 1);

            % Identity of the generative process
            if isnan(G{cond,sub}.Questions(2))
                text(S{sub}.Nstims/2, 3/4, 'Stochastic', 'HorizontalAlignment', ...
                    'Center', 'FontWeight', 'Bold');
            elseif G{cond,sub}.Questions(2) == 1
                text(S{sub}.Nstims/2, 3/4, 'Probabilistic', 'HorizontalAlignment', ...
                    'Center', 'FontWeight', 'Bold');
            elseif G{cond,sub}.Questions(2) == 2
                text(S{sub}.Nstims/2, 3/4, 'Deterministic', 'HorizontalAlignment', ...
                    'Center', 'FontWeight', 'Bold');
            end

            % Customize the axes
            axis([1,N,0,1]);
            set(gca, 'YTick', [], 'LineWidth', 2);

            % Display the generative process as a title
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            for i = [seq, seq + m, seq + m*2]
                subplot(3, m, i);
                if     reg == 3
                    title({'Fully', 'stochastic'});
                elseif reg == 1
                    title(sprintf('p(%s|%s) = %1.2f\np(%s|%s) = %1.2f', ...
                        letters{1}, letters{2}, G{cond,sub}.Rule(1), ...
                        letters{2}, letters{1}, G{cond,sub}.Rule(2)));
                elseif reg == 2
                    title({char(letters(G{cond,sub}.Rule))',''});
                end
            end
        end

        % Get the figure
        F{reg} = getframe(f);
        F{reg} = frame2im(F{reg});

        % Close the figure
        close(f);
    end

    % Save the figure
    fig = cell2mat(F');
    fig = fig(:,400:3550,:);
    filename = sprintf('F_Ind_Sub%02i.png', sub);
    pathtofile = fullfile(homedir, 'Finger tracking analyses', 'subfig', filename);
    imwrite(fig, pathtofile, 'png');
end
