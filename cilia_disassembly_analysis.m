clear all % Make sure to clear out everything before we start
%% 
% What this code does:
%% 
% # Categorizes cilia disassembly curve into three categories: Instant, Gradual, 
% and Combined (Gradual + Instant)
% # Determines the rate of disassembly (micron/minute as well as % decrease 
% per minute)
% # Normalizes curves for comparison
%% 
% How it works:
%% 
% # Smooth noise from raw data
% # Determines a "start point" for disassembly using a heuristic (where smoothed 
% derivative first becomes larger than the normal, slow change in length exhibited 
% by control curves, from a previous analysis)
% # Use start point to determine whether Gradual dissassembly occurred and the 
% disassembly rate
% # Use the drop in length at the very end to determine whether Instant disassembly 
% occurred.
%% 
% Terminology used in this code:
%% 
% * Terminal Deciliation (TD): period of contiguous cilia length reduction leading 
% to complete cilia loss
% * Instant Terminal Deciliation (ITD): When a cilium loses its entire length 
% between two consecutive measurements (and that loss is larger than a threshold 
% determined from control data). Possibly indicates a basal severing event. 
% * Gradual Terminal Deciliation (GTD): When a cilium's length reduces over 
% many consecutive measurements leading to loss. 
% * Both ITD and GTD can occur in the same observation, i.e., when GTD starts 
% before ITD occurs. This corresponds to the 'Combined' category. 

% Load data from an Excel file located in the same directory as this code.
[curves, labels] = xlsread("SampleDataset.xlsx", "Sheet 1");
labels = labels(1, :); % Just take the first row from labels, if necessary
%%
%%%% Helper variables %%%%

num_obs = round(size(curves, 2) / 2);   % # observations we are analyzing

% Initializers for result table
is_ITD = zeros(num_obs, 1);             % whether ITD occurred
is_GTD = zeros(num_obs, 1);             % whether GTD occurred
rate_ITD = zeros(num_obs, 1);           % rate of ITD (if present)
rate_GTD = zeros(num_obs, 1);           % rate of GTD (if present)
loss_pct_ITD = zeros(num_obs, 1);       % percent loss due to ITD
loss_pct_GTD = zeros(num_obs, 1);       % percent loss due to GTD
categorization = strings(num_obs, 2);   % category: I(TD), G(TD), C(ombined)
results = strings(num_obs, 11);         % result output table

% Thresholds from control analysis (NOTE: Make sure to input your own control data here!)
control_slope = -0.005;                 % on average control cilia reduce in length by this many mu/min
control_instant_stddev = 1.53;          % stddev of control data. we use this as a noise threshold 

% Initializers for normalized graphs
num_points = 1000;
normalized_instant = zeros(num_obs, num_points);
normalized_instant_labels = {};
normalized_gradual = zeros(num_obs, num_points);
normalized_gradual_labels = {};
normalized_combined = zeros(num_obs, num_points);
normalized_combined_labels = {};

columns = 1:2:size(curves,2);           % point to the "x" column of each cilium

for i = 1:numel(columns)
    %%% Data prep %%%
    col = columns(i);
    
    % get x,y data and clean it a bit
    
    obs = sortrows(curves(:, col:col+1)); % make sure they are sorted
    
    x = obs(:,1); % time  (min)
    y = obs(:,2); % length (um)
    
    x(isnan(x)) = []; % remove NaNs (import artifacts)
    y(isnan(y)) = []; % remove NaNs
    
    %%% Instant disassembly %%%
   
    % Simply check whether second-to-last y-value is larger than our threshold. 
    if (y(end-1) > control_instant_stddev)
        is_ITD(i) = 1;
    else
        is_ITD(i) = 0;
    end    
    
    %%% Gradual disassembly %%%
    
    % First, we need to smooth out the noisy data. We can provide the `smoothdata` function with a smoothing factor (default 0.25). 
    % We determined the best value for our data by testing a subset of our curves and comparing the result to an expert human opinion.
    smoothing_factor = 0.35;
    
    % We don't want instantaneous loss to skew our smoothing, 
    % so drop the last point if instant.
    if is_ITD(i) 
        maybe_drop_last = [x(1:end-1) y(1:end-1)];
        smooth_plot_x = x(1:end-1);
    else
        maybe_drop_last = [x y];
        smooth_plot_x = x;
    end
    
    % We had best results using a gaussian smoothing method
    [gauss_curve, smoothing_window] = smoothdata(maybe_drop_last, 'gaussian', 'SmoothingFactor', smoothing_factor);
    
    % Our dataset is a combination of both semi-automated Imaris measurements and manual measurements. These exhibit two distinct
    % noise levels which are accounted for with the following heuristic. 
    % If the window size is too small, barely any smoothing is actually done, so we bump the smoothing factor a bit and try again.
    
    % We chose this value for `min_window` based on our manual tests with a subset of our curves. 
    % It depends on the number of time points and the time interval.
    min_window = 6;
    if smoothing_window < min_window
        [gauss_curve, smoothing_window] = smoothdata(maybe_drop_last, 'gaussian', 'SmoothingFactor', smoothing_factor + 0.2);
    end
    
    smoothx = gauss_curve(:, 1);
    smoothy = gauss_curve(:, 2);
    
    % smooth derivatives. 
    dsmooth = diff(smoothy)./diff(smoothx);
    
    % check for descent using smoothed derivatives:
    %   get the last time index where derivative > control_slope
    %   any derivative after this point is negative, so the next index
    %   is the beginning of the last descent (if any).
    start_dsmooth_idx = find(dsmooth>control_slope, 1, 'last');
    
    % if we don't find any point that fits the above criteria, it means that the curve has been descending for 
    % the entire time, so we set `start_dsmooth` to beginning
    if isempty(start_dsmooth_idx)
        start_dsmooth = 1;
    else
        start_dsmooth = min(start_dsmooth_idx + 1, size(dsmooth, 1));
    end

    % is_gradual calc:
    %   if the start point we calculated above is before the last two points of the curve (i.e, not instant), 
    %   we say it's gradual
    
    if start_dsmooth < size(dsmooth, 1) - 1
        is_GTD(i) = 1;
    else
        is_GTD(i) = 0;
    end    
    
    %% categorization and category-related calculations
    
    % Normalization 
    interpolation_values = linspace(1, max(x), num_points);            % Normalize all time points to `num_points`
    interpolated_row = interp1(x, y, interpolation_values, 'linear');  % Interpolate values onto above
    normalized_row = interpolated_row./max(interpolated_row);          % Normalize length values such that max length in curve = 1
    
    % Here, we set this curve's label, rate and loss percent of ITD/GTD (if present) and populate normalized tables. 
    if is_ITD(i) && is_GTD(i)  % Combined
        categorization(i, :) = [string(labels{col}) "C"];
          
        rate_GTD(i) = (smoothy(end-1) - smoothy(start_dsmooth)) / (smooth_plot_x(end-1) - smooth_plot_x(start_dsmooth));
        rate_ITD(i) = smoothy(end)/(x(end) - x(end-1));
                
        loss_pct_GTD(i) = (smoothy(start_dsmooth) - smoothy(end-1)) / smoothy(start_dsmooth);
        loss_pct_ITD(i) = smoothy(end) / smoothy(start_dsmooth);
        
        normalized_combined(i, :) = normalized_row;
        normalized_combined_labels{end + 1} = labels{col};
    elseif is_ITD(i)
        categorization(i, :) = [string(labels{col}) "I"];
        
        rate_ITD(i) = smoothy(end)/(x(end) - x(end-1));
        loss_pct_ITD(i) = smoothy(end-1) / smoothy(start_dsmooth);
        
        normalized_instant(i, :) = normalized_row;
        normalized_instant_labels{end + 1} = labels{col};
    elseif is_GTD(i)
        categorization(i, :) = [string(labels{col}) "G"];
        
        rate_GTD(i) = (smoothy(end) - smoothy(start_dsmooth)) / (smooth_plot_x(end) - smooth_plot_x(start_dsmooth));
        loss_pct_GTD(i) = (smoothy(start_dsmooth) - smoothy(end)) / smoothy(start_dsmooth);
        
        normalized_gradual(i, :) = normalized_row;
        normalized_gradual_labels{end + 1} = labels{col};
    else
        categorization(i, :) = [string(labels{col}) "?"];
    end
        
    results(i, :) = [
        string(i),                                  % i
        string(labels{col}),                        % cilium #
        string(categorization(i, 2)),               % category
        string(x(start_dsmooth)),                   % TD start (min)
        string(x(end)),                             % end time (min)
        string(round(rate_GTD(i), 5)),              % GTD rate (um/min)
        string(round(rate_ITD(i), 5)),              % ITD rate (um/min)
        string(round(loss_pct_GTD(i), 5)*100),      % GTD pct loss
        string(round(loss_pct_ITD(i), 5)*100)       % ITD pct loss
        string(y(1)),                               % initial length (um)
        string(y(start_dsmooth))                    % length at TD start (um)
    ];

    % Plot a single graph with the following data:
    %    Raw data
    %    Smoothed data
    %    Start point from smooth data
    %    Start point projected onto raw data
    % Graphs are also labeled with calculated rates and lengths from above
    
    figure
    hold on
    plot(x, y, '-.o', 'DisplayName', labels{col})
    plot(smooth_plot_x, smoothy, '-x', 'DisplayName',['smooth ' mat2str(smoothing_window)])

    legend('Location','bestoutside');

    plot(smooth_plot_x(start_dsmooth), y(start_dsmooth), 'oc', 'MarkerSize',10, 'DisplayName', 'start (gauss)');
    h = plot(smooth_plot_x(start_dsmooth), smoothy(start_dsmooth), 'oc', 'MarkerSize',10, 'DisplayName', 'start (gauss)');
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    hold off
    % title
    result = join([
        'i = ' results(i, 1) ' | '...
        'cilium ' results(i, 2) ' | '...
        'category: ' results(i, 3) ' | '...
        'TD start ' results(i, 4) ' min | '...
        'end ' results(i, 5) ' min' newline ...
        'GTD rate ' results(i, 6) ' ?/min | '...
        'ITD rate ' results(i, 7) ' ?/min' newline...
        'GTD loss ' results(i, 8) '% | '...
        'ITD loss ' results(i, 9) '% | ' newline...
        'length (init) ' results(i, 10) '? | '...
        'length (TD start) ' results(i, 11) '? '
    ]);

    title(result)
    xlabel('Time (min)')
    ylabel('Length (um)')
end
% Create a table so we can write a CSV file with our results
resultsTable = array2table(results, 'VariableNames',{
    'i', 'cilia', 'category', 'td_start_time_min', 'end_time_min', ...
    'gtd_rate_um_per_min', 'itd_rate_um_per_min', ...
    'gtd_loss_pct', 'itd_loss_pct','length_i', 'length_td_start' 
})
writetable(resultsTable, 'results.csv')
%%
%% Plot normalized curves for each category

% remove all 0 rows from normalized data
normalized_combined = normalized_combined(any(normalized_combined,2),:);
normalized_instant = normalized_instant(any(normalized_instant,2),:);
normalized_gradual = normalized_gradual(any(normalized_gradual,2),:);

% generate mean curve

mean_combined = mean(normalized_combined);
mean_gradual = mean(normalized_gradual);
mean_instant = mean(normalized_instant);

% update labels
normalized_combined_labels{end + 1} = 'avg';
normalized_gradual_labels{end + 1} = 'avg';
normalized_instant_labels{end + 1} = 'avg';

% save normalized + average curves
combined_table = array2table([normalized_combined; mean_combined]', 'VariableNames', normalized_combined_labels)
gradual_table = array2table([normalized_gradual; mean_gradual]', 'VariableNames', normalized_gradual_labels)
instant_table = array2table([normalized_instant; mean_instant]', 'VariableNames', normalized_instant_labels)

writetable(combined_table, 'combined_normalized.csv')
writetable(instant_table, 'instant_normalized.csv')
writetable(gradual_table, 'gradual_normalized.csv')


% Plot the normalized curves
figure
plot(normalized_combined') % note the transposition
hold on
plot(mean_combined', 'LineWidth', 5)
hold off
title("Normalized Combined Data")
xlabel('Index')
ylabel('Normalized Length (um)')

figure
plot(normalized_gradual') % note the transposition
hold on
plot(mean_gradual', 'LineWidth',5)
hold off
title("Normalized Gradual Data")
xlabel('Index')
ylabel('Normalized Length (um)')

figure
plot(normalized_instant') % note the transposition
hold on
plot(mean_instant', 'LineWidth',5)
hold off
title("Normalized Instant Data")
xlabel('Index')
ylabel('Normalized Length (um)')