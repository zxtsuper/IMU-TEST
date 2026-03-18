% North Finding Demo Script
% Author: zxtsuper
% Date: 2026-03-18
% This script demonstrates the process of finding north based on IMU data.
% The primary output is the y-axis heading from North (CW+).

% Load the IMU data
load('imu_data.mat');

% Process the data
res = process_imu_data(imu_data);

% Calculate y-axis heading (CW+)
heading_y_deg = mod(res.heading_deg + 90, 360);

% Plot the results
figure;
plot(res.time, heading_y_deg, 'LineWidth', 2);
hold on;

% Calculate 1-sigma band
y_sigma = std(res.yaw_history_deg);
plot(res.time, mod(res.yaw_history_deg + 90, 360) + y_sigma, 'r--');
plot(res.time, mod(res.yaw_history_deg + 90, 360) - y_sigma, 'r--');

% Figure attributes
title('Y-Axis Heading from North (CW+)');
xlabel('Time (s)');
ylabel('Heading (degrees)');
grid on;
legend('Y-Axis Heading', '1-Sigma Band');
hold off;

% Save figure
saveas(gcf, 'y_axis_heading_plot.png');
