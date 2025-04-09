clc; 
clear; 
close all; 

cmp = lines(8);

%% Data Reading and Parsing

residue_pairs = load('sb_probabilities_simple.dat');

%% Residue Groups
group1_range = 33:44;   % Range of the first group of residues
group2_range = 90:422;  % Range of the second group of residues

%% Sum of Probability Forming Salt Bridges
% 1st
[group1_res, ~, idx1] = unique(residue_pairs(:,1));
group1_probs = accumarray(idx1, residue_pairs(:,3), [], @sum);

% 2nd
[group2_res, ~, idx2] = unique(residue_pairs(:,2));
group2_probs = accumarray(idx2, residue_pairs(:,3), [], @sum);

%% Data Sorting and Alignment
% Ensure display of all residues within the defined ranges
full_group1 = group1_range(:);
[~, loc1] = ismember(full_group1, group1_res);
group1_final = [full_group1, zeros(length(full_group1),1)];
for i = 1:length(full_group1)
    if loc1(i) ~= 0
        group1_final(i,2) = group1_probs(loc1(i));
    end
end

full_group2 = group2_range(:);
[~, loc2] = ismember(full_group2, group2_res);
group2_final = [full_group2, zeros(length(full_group2),1)];
for i = 1:length(full_group2)
    if loc2(i) ~= 0
        group2_final(i,2) = group2_probs(loc2(i));
    end
end

%% Output
% 1st
fid = fopen('group1_summary.csv', 'w');
fprintf(fid, 'Residue,Total_Probability\n');
fprintf(fid, '%d,%.2f\n', group1_final');
fclose(fid);

% 2nd
fid = fopen('group2_summary.csv', 'w');
fprintf(fid, 'Residue,Total_Probability\n');
fprintf(fid, '%d,%.2f\n', group2_final');
fclose(fid);

%% Visualization
figure('Position', [100, 100, 1200, 500], 'Color', 'w')

% 1st Group
subplot(1,2,1)
bar(group1_final(:,1), group1_final(:,2), 'FaceColor', [0.2 0.6 0.8])
title('Group1 (33-44) Residues Salt Bridge Probability', 'FontSize',12, 'FontWeight','bold')
xlabel('Residue Number', 'FontSize',10)
ylabel('Total Probability (%)', 'FontSize',10)
xlim([32 45])
grid on
set(gca, 'XTick', group1_range, 'XTickLabelRotation',90, 'FontSize',8)

% 2nd Group
subplot(1,2,2)
bar(group2_final(:,1), group2_final(:,2), 'FaceColor', [0.8 0.4 0.1])
title('Group2 (90-422) Residues Salt Bridge Probability', 'FontSize',12, 'FontWeight','bold')
xlabel('Residue Number', 'FontSize',10)
xlim([40 423])
grid on
set(gca, 'XTick', 90:50:420, 'XTickLabelRotation',45, 'FontSize',8)

% Save the figure
exportgraphics(gcf, 'saltbridge_summary.png', 'Resolution',300)

%% Combine Data with Probability > 10
group1_data = readtable('group1_summary.csv');
group2_data = readtable('group2_summary.csv');

group1_data.Group = repmat("Group1", height(group1_data), 1);
group2_data.Group = repmat("Group2", height(group2_data), 1);
combined_data = [group1_data; group2_data];

high_prob = combined_data(combined_data.Total_Probability > 5, :);

writetable(high_prob, 'sb.dat', ...
          'Delimiter', 'tab', ...
          'WriteVariableNames', true, ...
          'FileType', 'text');