% NORTH_FINDING_DEMO  MEMS-IMU gyrocompass north-finding demo (MATLAB R2020).
%
% Purpose
% -------
%   Demonstrates how to use IMU data alone (accelerometers + gyroscopes)
%   to estimate the compass heading (yaw / north-finding) with an EKF.
%   No GNSS, magnetometer, or vision sensors are used.
%
% Hardware assumptions
% --------------------
%   - MEMS IMU with gyro bias / drift grade ~5 deg/h
%   - Static (stationary) placement during measurement
%   - Body z-axis points roughly DOWN (acc_z ≈ -g at rest)
%
% Coordinate convention
% ----------------------
%   Navigation frame : ENU  (x = East, y = North, z = Up)
%   Rotation Cbn = Rx(phi) * Ry(theta) * Rz(psi)  maps ENU → body
%   Heading output : atan2d(cos(psi), -sin(psi))   [deg, North = 0, CW+]
%   This is the compass bearing of the IMU x-axis (column 2 of data).
%
% Usage
% -----
%   1. Place the IMU file at  ../data/test.txt  (or update DATA_FILE below).
%   2. Set your geodetic latitude (LATITUDE_DEG).
%   3. Run this script in MATLAB R2020 or later.
%
% Dependencies
% ------------
%   No extra toolboxes required for core functionality.
%   If Signal Processing Toolbox is present, the low-pass filter uses
%   butter + filtfilt (better roll-off); otherwise a moving-average is used.
%
% Output
% ------
%   Console: heading estimate with 1-sigma confidence.
%   Figures: (1) raw sensor data, (2) static detection, (3) EKF convergence.

clear; close all; clc;

% ==========================================================================
%  USER CONFIGURATION — edit these values
% ==========================================================================

DATA_FILE     = '../data/test.txt';  % path relative to this script
LATITUDE_DEG  = 31.23;               % your geodetic latitude [degrees N]
                                     % e.g. 39.9 = Beijing, 31.2 = Shanghai

% Gyro / acc unit in the data file
GYRO_UNIT     = 'rad_s';   % 'rad_s' or 'deg_s'
ACC_UNIT      = 'm_s2';    % 'm_s2'  or 'g'

% EKF noise parameters  (tune to your sensor)
SIGMA_GYRO    = 2.0e-5;    % post-lowpass gyro noise STD [rad/s]
                            % Approx: raw_gyro_noise / sqrt(fs/lp_cutoff)
                            % e.g. 2e-4 / sqrt(200/1.5) ≈ 2e-5 rad/s
SIGMA_BIAS_RW = 2.0e-6;    % gyro bias random-walk [rad/s/sqrt(s)]
SIGMA_YAW_RW  = 1.0e-5;    % yaw random-walk [rad/sqrt(s)]
MAX_BIAS      = 3.0e-4;    % expected max bias magnitude [rad/s]
                            % 5 deg/h ≈ 2.42e-4 rad/s; set slightly above

% Low-pass filter
LP_CUTOFF_HZ  = 1.5;        % 0 = disable; 0.5–3 Hz is typical for static use

% Window length limits
MIN_STATIC_S  = 10;          % minimum static segment length [s]
MAX_STATIC_S  = 60;          % maximum; cap at 60–120 s to limit temp drift
MAX_EKF_S     = 60;          % maximum data fed to EKF [s]

% ==========================================================================
%  STEP 1: Load data
% ==========================================================================
script_dir = fileparts(mfilename('fullpath'));
data_path  = fullfile(script_dir, DATA_FILE);

fprintf('=== North-Finding Demo ===\n');
fprintf('Loading IMU data from: %s\n', data_path);

read_opts.gyro_unit = GYRO_UNIT;
read_opts.acc_unit  = ACC_UNIT;
read_opts.verbose   = true;

imu = read_imu_csv(data_path, read_opts);

t    = imu.t;
acc  = imu.acc;
gyro = imu.gyro;
fs   = imu.fs;
N    = imu.N;

% ==========================================================================
%  STEP 2: Static segment detection
% ==========================================================================
fprintf('\n--- Static Segment Detection ---\n');

static_opts.min_static_s  = MIN_STATIC_S;
static_opts.max_static_s  = MAX_STATIC_S;
static_opts.verbose       = true;

[i0, i1, stat_score] = detect_static_segments(t, acc, gyro, static_opts);

t_static   = t(i0:i1);
acc_static = acc(i0:i1, :);
gyr_static = gyro(i0:i1, :);

fprintf('Static window duration: %.2f s  (samples %d–%d)\n', ...
    t_static(end) - t_static(1), i0, i1);

% ==========================================================================
%  STEP 3: Coarse alignment — roll and pitch from gravity
% ==========================================================================
fprintf('\n--- Coarse Alignment (roll / pitch from gravity) ---\n');

g0       = 9.80665;
acc_mean = mean(acc_static, 1)';
acc_mag  = norm(acc_mean);

u     = acc_mean / max(acc_mag, 0.5 * g0);
roll  = atan2(-u(2), u(3));
pitch = atan2( u(1), sqrt(u(2)^2 + u(3)^2));

fprintf('Mean |acc| = %.4f m/s^2  (g0 = %.4f)\n', acc_mag, g0);
fprintf('Roll  = %.3f deg\n', rad2deg(roll));
fprintf('Pitch = %.3f deg\n', rad2deg(pitch));

if abs(acc_mag - g0) > 0.5
    warning(['|acc| differs from g0 by %.2f m/s^2. ' ...
             'Verify acc_unit and that the device is stationary.'], ...
        abs(acc_mag - g0));
end

% ==========================================================================
%  STEP 4: North-finding / heading estimation (EKF)
% ==========================================================================
fprintf('\n--- EKF North-Finding ---\n');

ekf_opts.max_seconds     = MAX_EKF_S;
ekf_opts.lp_cutoff_hz    = LP_CUTOFF_HZ;
ekf_opts.sigma_gyro_meas = SIGMA_GYRO;
ekf_opts.sigma_bias_rw   = SIGMA_BIAS_RW;
ekf_opts.sigma_yaw_rw    = SIGMA_YAW_RW;
ekf_opts.max_bias_rad_s  = MAX_BIAS;
ekf_opts.att_update_hz   = 2;
ekf_opts.verbose         = true;

res = estimate_yaw_northfinding( ...
    t_static, acc_static, gyr_static, LATITUDE_DEG, ekf_opts);

% ==========================================================================
%  STEP 5: Print summary
% ==========================================================================
fprintf('\n=== RESULT ===\n');
fprintf('Heading (body x-axis from North, CW+) = %.2f deg\n', res.heading_deg);
fprintf('1-sigma uncertainty                   = %.2f deg\n', ...
    sqrt(max(res.cov_psi_deg2, 0)));
fprintf('Roll  = %.3f deg\n', res.roll_deg);
fprintf('Pitch = %.3f deg\n', res.pitch_deg);
fprintf('Gyro bias estimate [rad/s]: gx=%.3e  gy=%.3e  gz=%.3e\n', ...
    res.bias_est_rad_s(1), res.bias_est_rad_s(2), res.bias_est_rad_s(3));
fprintf('Duration processed: %.2f s\n', res.used_duration_s);
fprintf('\nNote: For a 5 deg/h MEMS gyro, expect ~5–20 deg accuracy after\n');
fprintf('30–60 s of quiet static data. Longer windows improve SNR.\n');

% ==========================================================================
%  STEP 6: Visualisation
% ==========================================================================
fprintf('\n--- Plotting ---\n');
t_rel  = t - t(1);
ts_rel = t_static - t_static(1);

% ---- Figure 1: Raw sensor data ----------------------------------------
figure('Name', 'Raw IMU data', 'NumberTitle', 'off');

subplot(2,1,1);
plot(t_rel, acc(:,1), 'r', t_rel, acc(:,2), 'g', t_rel, acc(:,3), 'b');
hold on;
h_patch = patch([t(i0)-t(1), t(i1)-t(1), t(i1)-t(1), t(i0)-t(1)], ...
                [min(acc(:))-1, min(acc(:))-1, max(acc(:))+1, max(acc(:))+1], ...
                [0.8 0.95 0.8], 'FaceAlpha', 0.4, 'EdgeColor', 'none'); %#ok<NASGU>
grid on; ylabel('Acc (m/s^2)');
legend('ax','ay','az','Location','best');
title('Accelerometer (green patch = static window)');

subplot(2,1,2);
plot(t_rel, gyro(:,1)*1e3, 'r', t_rel, gyro(:,2)*1e3, 'g', t_rel, gyro(:,3)*1e3, 'b');
hold on;
ylims = ylim;
patch([t(i0)-t(1), t(i1)-t(1), t(i1)-t(1), t(i0)-t(1)], ...
      [ylims(1), ylims(1), ylims(2), ylims(2)], ...
      [0.8 0.95 0.8], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
grid on; ylabel('Gyro (\times10^{-3} rad/s)');
legend('gx','gy','gz','Location','best');
title('Gyroscope');
xlabel('Time (s)');

% ---- Figure 2: Sensor magnitudes and static detection -----------------
figure('Name', 'Static detection', 'NumberTitle', 'off');

acc_norm  = sqrt(sum(acc.^2, 2));
gyro_norm = sqrt(sum(gyro.^2, 2));

subplot(2,1,1);
plot(t_rel, acc_norm, 'b');
hold on;
plot([t(i0)-t(1), t(i0)-t(1)], ylim, 'r--', ...
     [t(i1)-t(1), t(i1)-t(1)], ylim, 'r--');
yline(g0, 'k:', 'g_0');
grid on; ylabel('|acc| (m/s^2)');
legend('|acc|','window bounds','Location','best');
title(sprintf('Accelerometer magnitude (static score = %.4g)', stat_score));

subplot(2,1,2);
plot(t_rel, gyro_norm * 1e4, 'm');
hold on;
plot([t(i0)-t(1), t(i0)-t(1)], ylim, 'r--', ...
     [t(i1)-t(1), t(i1)-t(1)], ylim, 'r--');
grid on; ylabel('|gyro| (\times10^{-4} rad/s)');
legend('|gyro|','window bounds','Location','best');
title('Gyroscope magnitude');
xlabel('Time (s)');

% ---- Figure 3: EKF heading convergence --------------------------------
figure('Name', 'EKF Heading Convergence', 'NumberTitle', 'off');

M      = numel(res.yaw_history_deg);
t_ekf  = ts_rel(res.used_indices);

subplot(2,1,1);
plot(t_ekf, res.yaw_history_deg, 'b', 'LineWidth', 1.5);
hold on;
sigma_hist = sqrt(max(res.cov_history_deg2, 0));
fill([t_ekf; flipud(t_ekf)], ...
     [res.yaw_history_deg + sigma_hist; flipud(res.yaw_history_deg - sigma_hist)], ...
     'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
grid on; ylabel('Heading (deg)');
title(sprintf('EKF Heading Estimate  (final = %.2f ± %.2f deg)', ...
    res.heading_deg, sqrt(max(res.cov_psi_deg2,0))));
legend('heading','1-\sigma bound','Location','best');

subplot(2,1,2);
plot(t_ekf, res.bias_history * 1e4, 'LineWidth', 1.2);
grid on; ylabel('Gyro bias (\times10^{-4} rad/s)');
legend('b_x','b_y','b_z','Location','best');
title('Estimated gyro bias (EKF)');
xlabel('Time (s)');

% ---- Figure 4: EKF innovation -----------------------------------------
figure('Name', 'EKF Innovation', 'NumberTitle', 'off');
plot(t_ekf, res.innovation * 1e4, 'LineWidth', 1);
grid on;
ylabel('Innovation (\times10^{-4} rad/s)');
xlabel('Time (s)');
legend('\nu_x','\nu_y','\nu_z','Location','best');
title('EKF Innovation Sequence (should approach 0 as EKF converges)');

fprintf('Done. Figures plotted.\n');
