% Update the default console output and convergence plot focus on heading of the body y-axis.

heading_y_deg = mod(res.heading_deg + 90, 360);
fprintf('Primary Heading (Y-axis): %f degrees\n', heading_y_deg);

% Update figure 3 to plot heading_y history and its 1-sigma bound
plot(mod(yaw_history_deg + 90, 360));
hold on;
% Code for 1-sigma bounds not shown in prompt
% ... 

% Finish plotting code here
hold off;