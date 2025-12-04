clc; clear; close all; box on;
hold on;

% --- 1. Data Path Settings ---
% Paths for six parallel simulation replicas (run1 to run6)
WT={'/run1_water_order_PATH/';...
    '/run2_water_order_PATH/';...
    '/run3_water_order_PATH/';...
    '/run4_water_order_PATH/';...
    '/run5_water_order_PATH/';...
    '/run6_water_order_PATH/'};

filename = 'wat-in-pore-order.dat'; 

% --- 2. Data Reading and Concatenation ---
systemData_or = []; % Initialize array to store all Order Parameters (P1)
systemData_z = [];  % Initialize array to store all Z-coordinates

for j = 1:length(WT)
    filePath = strcat(WT{j}, filename); 

    data = importdata(filePath); 
    
    % Logical indexing: filter data where the first column is > 0 
    % (usually filtering out invalid frames or negative time steps)
    logical_idx = data(:,1) > 0;
    t = data(logical_idx,1);            % Time/Frame
    z = data(logical_idx, 2);           % Z-coordinate
    order_param = data(logical_idx, 3); % Corresponding Order Parameter (P1)
    
    % Concatenate data from the current run into the aggregate arrays
    systemData_or = [systemData_or; order_param];
    systemData_z = [systemData_z; z];
end

% --- 3. Data Binning and Statistics ---
num_bins = 40;  
z_edges = linspace(-20.5, 20.5, 42); % Define Z-axis bin edges, covering range -20.5 to 20.5 Å

% Map all Z-coordinates to bins; bin_idx stores the bin index for each data point
[~, ~, bin_idx] = histcounts(systemData_z, z_edges);

% Calculate the mean order parameter for each bin
% 'accumarray' groups systemData_or by bin_idx and applies the mean function
bin_means = accumarray(bin_idx(bin_idx>0), systemData_or(bin_idx>0), [length(z_edges)-1, 1], @mean);

% Define bin centers for plotting (step size = 1 Å)
bin_centers = [-20:1:20];

% Calculate standard deviation (std) for each bin (optional, for reference)
bin_std = accumarray(bin_idx(bin_idx>0), systemData_or(bin_idx>0), [length(z_edges)-1, 1], @std);

% --- 4. Construct Boxchart Data Matrix ---
% 'boxchart' requires each column to represent a group. 
% Since the number of data points per bin varies, we must construct a NaN-padded matrix.

% 4.1 Find the maximum number of data points in a single bin to determine matrix row count
max_length = 0;
for bin = 1:(length(z_edges)-1)
    current_length = sum(bin_idx == bin);
    if current_length > max_length
        max_length = current_length;
    end
end

% 4.2 Initialize matrix (rows=max_length, cols=num_bins), pre-filled with NaN
bin_data_matrix = NaN(max_length, length(z_edges)-1);

% 4.3 Fill the matrix columns with actual data from each bin
for bin = 1:(length(z_edges)-1)
    current_data = systemData_or(bin_idx == bin);
    if ~isempty(current_data)
        bin_data_matrix(1:numel(current_data), bin) = current_data;
    end
end

% --- 5. Plotting ---
cmp = lines(1); % Use MATLAB's default color scheme (first color)
a = bin_centers; 

% Draw boxchart
% 'WhiskerLineStyle' set to 'none' to hide whiskers (showing only IQR)
boxchart(bin_data_matrix,'LineWidth',1.8,'WhiskerLineStyle', 'none')

% Overlay mean value curve
p = plot(bin_means, '-o');
p.LineWidth = 2;
p.MarkerSize = 10;

% Axes and Style Adjustments ---
xlabel('z (Å)');
ylabel('Order parameter P1');
set(gca, 'YTick', -1:0.5:1) % Set Y-axis ticks
ylim([-1 1]) % Fix Y-axis limits

% Custom X-axis labels: Prevent label overcrowding
labels = a;
labels_str = strings(size(labels));
% Logic: Display label every 5 points, ensure the last point is also shown
labels_str([1:5:end,end]) = string(labels([1:5:end,end]));
set(gca, 'XTickLabel', labels_str, 'XTickLabelRotation', 0);

% Grid settings
ax = gca;  
ax.YGrid = 'on';        
ax.YMinorGrid = 'on';  
ax.GridLineWidth = 1.5;  
ax.GridAlpha = 0.2;  

% Font Size Control ---
set(gca, 'linewidth', 3, 'fontweight', 'normal', 'fontname', 'Arial');

ax.XAxis.FontSize = 35;  
ax.YAxis.FontSize = 35;  

xlabel('z (Å)', 'FontSize', 40);              
ylabel('Order parameter P1', 'FontSize', 40); 

% Image Output ---
set(gcf,'position',[100,100,1600,520]); 
saveas(gcf,'p1_boxchart.png') 