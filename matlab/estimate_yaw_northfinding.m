function result = estimate_yaw_northfinding(t, acc, gyro, lat_deg, opts)
% ESTIMATE_YAW_NORTHFINDING  EKF-based gyrocompass / north-finding.
%
%   RESULT = ESTIMATE_YAW_NORTHFINDING(T, ACC, GYRO, LAT_DEG)
%   RESULT = ESTIMATE_YAW_NORTHFINDING(T, ACC, GYRO, LAT_DEG, OPTS)
%
%   Uses a 4-state Extended Kalman Filter to estimate the heading (yaw)
%   and gyroscope constant bias from static IMU data.  The Earth's
%   rotation rate is used as the observable "signal" to determine North.
%
% -------------------------------------------------------------------------
%   COORDINATE CONVENTION
% -------------------------------------------------------------------------
%   Navigation frame : ENU  (East = x, North = y, Up = z)
%   Body frame       : generic; body-z should point roughly down for
%                      the default roll/pitch estimation to work (i.e.
%                      acc_z < 0 when the sensor sits on a flat table).
%
%   Rotation (ENU -> body):
%       Cbn = Rx(phi) * Ry(theta) * Rz(psi)
%   where psi (yaw) is the EKF state, phi/theta (roll/pitch) are
%   derived from the accelerometer gravity direction.
%
%   Physical meaning of yaw output:
%       heading_deg = atan2d(cos(psi), -sin(psi))  [degrees from North, CW+]
%   This is the compass bearing of the body x-axis.
%
% -------------------------------------------------------------------------
%   EKF STATE  x = [psi; b_gx; b_gy; b_gz]  (4 x 1)
% -------------------------------------------------------------------------
%   Process model (static / random walk):
%       psi_{k+1}  = psi_k  + w_psi
%       b_g_{k+1}  = b_g_k  + w_b
%
%   Measurement model:
%       z_k = gyro_body_k
%       h(x) = Cbn(psi, roll, pitch) * w_ie_ENU + b_g
%       H    = [d(Cbn*w_ie)/d(psi), I_3x3]
%
% -------------------------------------------------------------------------
%   INPUTS
% -------------------------------------------------------------------------
%   T       - Nx1 timestamp vector (seconds)
%   ACC     - Nx3 accelerometer matrix [ax, ay, az] (m/s^2)
%   GYRO    - Nx3 gyroscope matrix    [gx, gy, gz] (rad/s)
%   LAT_DEG - geodetic latitude (degrees North); default 31.23 deg
%   OPTS    - optional configuration struct (see defaults below)
%
% -------------------------------------------------------------------------
%   OPTS FIELDS  (all optional, with sensible defaults)
% -------------------------------------------------------------------------
%   max_seconds      - maximum window length to process (default 180 s)
%   lp_cutoff_hz     - low-pass filter cutoff (default 1.5 Hz)
%                      Set 0 to disable; uses simple moving average
%                      (no toolbox required).
%   sigma_gyro_meas  - gyro measurement noise std (default 2e-4 rad/s)
%                      Tune to actual gyro noise after low-pass filtering.
%   sigma_bias_rw    - gyro bias random-walk std (default 2e-6 rad/s/sqrt(s))
%   sigma_yaw_rw     - yaw random-walk std (default 1e-5 rad/sqrt(s))
%                      Small value; keeps yaw from drifting without cause.
%   att_update_hz    - rate at which roll/pitch are re-estimated from
%                      accelerometers (default 2 Hz)
%   trim_percent     - percentage of samples trimmed at each end for the
%                      robust mean used in roll/pitch init (default 20)
%   bias_init        - 3x1 initial gyro bias estimate (default [0;0;0])
%   yaw_init         - initial yaw guess in radians (default 0)
%   verbose          - print progress (default true)
%
% -------------------------------------------------------------------------
%   OUTPUTS  (struct fields)
% -------------------------------------------------------------------------
%   yaw_deg          - final heading estimate, degrees from North CW+
%   heading_deg      - same as yaw_deg (alias for clarity)
%   yaw_history_deg  - Mx1 heading history over the processing window
%   roll_deg         - final roll  estimate (degrees)
%   pitch_deg        - final pitch estimate (degrees)
%   bias_est_rad_s   - 3x1 estimated gyro bias (rad/s)
%   bias_history     - Mx3 bias history (rad/s)
%   innovation       - Mx3 EKF innovation sequence
%   cov_psi_deg2     - final yaw covariance (deg^2); sqrt = 1-sigma bound
%   used_duration_s  - duration of data processed (seconds)
%   used_indices     - indices into the input arrays that were processed
%   note             - text summary of the run
%
% -------------------------------------------------------------------------
%   ACCURACY GUIDANCE (5 deg/h MEMS gyro)
% -------------------------------------------------------------------------
%   Gyro bias dominates; typical 1-sigma heading accuracy:
%     - After  30 s:  ~15-30 deg (bias barely visible above noise)
%     - After  60 s:  ~ 5-15 deg
%     - After 120 s:  ~ 3- 8 deg
%   Accuracy improves with lower noise floor; 180 s is the recommended cap
%   to limit temperature-drift effects.
%
% -------------------------------------------------------------------------
if nargin < 4 || isempty(lat_deg)
    lat_deg = 31.23;
    warning('estimate_yaw_northfinding: latitude not specified; using default %.2f deg N.', lat_deg);
end
if nargin < 5 || isempty(opts)
    opts = struct();
end

% --- Defaults ---
opts = yn_setdef(opts, 'max_seconds',     180);
opts = yn_setdef(opts, 'lp_cutoff_hz',    1.5);
% Post-lowpass measurement noise: at 200 Hz with 1.5 Hz cut-off the box
% filter reduces white noise by ~sqrt(fs/fc) ≈ 11x.
% If your raw gyro noise is ~2e-4 rad/s, after filtering sigma ≈ 2e-5 rad/s.
% Tune this to std(gyro_filtered) from a quiet static segment.
opts = yn_setdef(opts, 'sigma_gyro_meas', 2.0e-5);
opts = yn_setdef(opts, 'sigma_bias_rw',   2.0e-6);
opts = yn_setdef(opts, 'sigma_yaw_rw',    1.0e-5);
% max_bias_rad_s: expected maximum gyro bias magnitude (1-sigma prior).
% For a 5 deg/h IMU: 5*(pi/180)/3600 ≈ 2.4e-4 rad/s.  Set slightly
% above worst-case to avoid over-constraining the EKF.
opts = yn_setdef(opts, 'max_bias_rad_s',  3.0e-4);
opts = yn_setdef(opts, 'att_update_hz',   2);
opts = yn_setdef(opts, 'trim_percent',    20);
opts = yn_setdef(opts, 'bias_init',       zeros(3,1));
opts = yn_setdef(opts, 'yaw_init',        0);
opts = yn_setdef(opts, 'verbose',         true);

% ==========================================================================
%  1. Preprocessing
% ==========================================================================
t    = t(:);
N    = numel(t);
dt0  = median(diff(t));
fs   = 1 / dt0;

% Limit to max_seconds
i_end = find(t <= t(1) + opts.max_seconds, 1, 'last');
if isempty(i_end)
    i_end = N;
end
idx  = (1:i_end)';

tw   = t(idx);
accw = acc(idx, :);
grw  = gyro(idx, :);

% Low-pass filter (noise suppression) — no toolbox required
accw_f = yn_lowpass(accw, fs, opts.lp_cutoff_hz);
grw_f  = yn_lowpass(grw,  fs, opts.lp_cutoff_hz);

% ==========================================================================
%  2. Earth rotation rate in ENU
% ==========================================================================
w_ie_n = earth_rate(lat_deg, 'ENU');   % [0; Om*cos(lat); Om*sin(lat)]

% ==========================================================================
%  3. Initial roll / pitch from accelerometer
% ==========================================================================
init_len = min(round(2 * fs), size(accw_f, 1));
f0       = yn_robust_mean(accw_f(1:init_len, :), opts.trim_percent);
[roll0, pitch0] = yn_roll_pitch(f0);

% ==========================================================================
%  4. EKF initialisation
% ==========================================================================
% State: x = [psi; b_gx; b_gy; b_gz]
x = [opts.yaw_init; opts.bias_init(:)];

sg = opts.sigma_gyro_meas;
R  = (sg^2) * eye(3);

sb = opts.sigma_bias_rw;
sy = opts.sigma_yaw_rw;
Q  = diag([ sy^2, sb^2, sb^2, sb^2 ]);   % process noise PSD

% Initial covariance: yaw unknown, bias constrained by sensor spec.
% A tighter bias prior (based on sensor grade) is CRITICAL for
% observability: the bias prior must be comparable to the earth-rate
% signal (~6e-5 rad/s at mid-latitudes).
mb = opts.max_bias_rad_s;
P  = diag([ (pi)^2, mb^2, mb^2, mb^2 ]);

% ==========================================================================
%  5. EKF main loop
% ==========================================================================
M        = numel(tw);
yaw_hist = zeros(M, 1);
bg_hist  = zeros(M, 3);
innov    = zeros(M, 3);
cov_hist = zeros(M, 1);

att_step = max(1, round(fs / opts.att_update_hz));
roll     = roll0;
pitch    = pitch0;

if opts.verbose
    fprintf('[northfinding] Processing %.1f s of data at %.1f Hz\n', ...
        tw(end)-tw(1), fs);
    fprintf('[northfinding] Init roll=%.2f deg  pitch=%.2f deg\n', ...
        rad2deg(roll), rad2deg(pitch));
end

for k = 1:M
    % ---- time step
    if k > 1
        dtk = tw(k) - tw(k-1);
        if dtk <= 0 || dtk > 1
            dtk = dt0;   % guard against bad timestamps
        end
    else
        dtk = dt0;
    end

    % ---- predict (random walk process model, F=I)
    Qd = Q * dtk;
    P  = P + Qd;

    % ---- refresh roll/pitch from accelerometer periodically
    if mod(k - 1, att_step) == 0
        wlen  = max(5, round(0.5 * fs));
        ks    = max(1, k - wlen + 1);
        fbar  = yn_robust_mean(accw_f(ks:k, :), opts.trim_percent);
        [roll, pitch] = yn_roll_pitch(fbar);
    end

    % ---- measurement update
    z   = grw_f(k, :)';          % 3x1 gyro measurement
    psi = x(1);
    bg  = x(2:4);

    Cbn  = yn_Cbn(psi, pitch, roll);
    h    = Cbn * w_ie_n + bg;    % predicted gyro reading (3x1)

    % Jacobian of h w.r.t. psi (numerical differentiation)
    eps2     = 1e-7;
    Cbn2     = yn_Cbn(psi + eps2, pitch, roll);
    dh_dpsi  = (Cbn2 * w_ie_n - Cbn * w_ie_n) / eps2;   % 3x1

    H = [dh_dpsi, eye(3)];    % 3x4

    inn = z - h;              % innovation
    S   = H * P * H' + R;
    K   = P * H' / S;

    x = x + K * inn;
    % Joseph stabilised covariance update (numerically robust)
    ImKH = eye(4) - K * H;
    P    = ImKH * P * ImKH' + K * R * K';

    % wrap yaw to [-pi, pi]
    x(1) = yn_wrap_pi(x(1));

    yaw_hist(k)  = x(1);
    bg_hist(k,:) = x(2:4)';
    innov(k,:)   = inn';
    cov_hist(k)  = P(1,1);
end

% ==========================================================================
%  6. Output
% ==========================================================================
psi_final    = x(1);
heading_deg  = yn_psi_to_heading(psi_final);

result.yaw_deg          = heading_deg;
result.heading_deg      = heading_deg;
result.yaw_history_deg  = yn_psi_to_heading_vec(yaw_hist);
result.roll_deg         = rad2deg(roll);
result.pitch_deg        = rad2deg(pitch);
result.bias_est_rad_s   = x(2:4);
result.bias_history     = bg_hist;
result.innovation       = innov;
result.cov_psi_deg2     = cov_hist(end) * (180/pi)^2;
result.cov_history_deg2 = cov_hist * (180/pi)^2;
result.used_duration_s  = tw(end) - tw(1);
result.used_indices     = idx;
result.note = sprintf( ...
    ['EKF gyrocompass | lat=%.2f N | %.1f s processed | ' ...
     'heading=%.2f deg (1-sigma ~%.2f deg)'], ...
    lat_deg, result.used_duration_s, heading_deg, ...
    sqrt(max(result.cov_psi_deg2, 0)));

if opts.verbose
    fprintf('[northfinding] Heading = %.2f deg from North (1-sigma ~%.2f deg)\n', ...
        heading_deg, sqrt(max(result.cov_psi_deg2, 0)));
    fprintf('[northfinding] Gyro bias est (rad/s): [%.3e, %.3e, %.3e]\n', ...
        x(2), x(3), x(4));
end
end

% ==========================================================================
%  Local helpers (all no-toolbox)
% ==========================================================================

function Cbn = yn_Cbn(psi, theta, phi)
% Rotation from ENU nav to body: Cbn = Rx(phi)*Ry(theta)*Rz(psi)
cp = cos(psi);  sp = sin(psi);
ct = cos(theta);st = sin(theta);
cf = cos(phi);  sf = sin(phi);

Rz = [cp, -sp, 0;  sp,  cp, 0;  0, 0, 1];
Ry = [ct,  0, st;   0,   1, 0; -st, 0, ct];
Rx = [ 1,  0,  0;   0,  cf,-sf;  0, sf, cf];
Cbn = Rx * Ry * Rz;
end

function [roll, pitch] = yn_roll_pitch(f_body)
% Roll and pitch from specific force vector in body frame.
% f_body = Cbn * [0;0;g0]  =>  f_body / g0 = [sin(theta);
%                                               -sin(phi)*cos(theta);
%                                                cos(phi)*cos(theta)]
g0 = 9.80665;
u  = f_body(:) / max(norm(f_body), g0 * 0.5);   % guard against near-zero
pitch = atan2(u(1), sqrt(u(2)^2 + u(3)^2));
roll  = atan2(-u(2), u(3));
end

function heading = yn_psi_to_heading(psi)
% Convert EKF yaw state psi to compass heading [deg, North=0, CW+].
% body-x in ENU = [cos(psi); -sin(psi); 0]  (derived from Cbn convention)
% heading = atan2(East_component, North_component)
%         = atan2(cos(psi), -sin(psi))
heading = mod(atan2d(cos(psi), -sin(psi)), 360);
end

function h_vec = yn_psi_to_heading_vec(psi_vec)
% Vectorised version of yn_psi_to_heading.
h_vec = mod(atan2(cos(psi_vec), -sin(psi_vec)) * (180/pi), 360);
end

function psi = yn_wrap_pi(psi)
% Wrap angle to [-pi, pi].
psi = mod(psi + pi, 2*pi) - pi;
end

function Xf = yn_lowpass(X, fs, fc)
% Zero-phase low-pass filter using a simple moving-average (no toolbox).
% fc: cutoff frequency in Hz; the window length = round(fs/fc) samples.
% Falls back to Signal Processing Toolbox filtfilt/butter if available.
if fc <= 0 || fc >= fs / 2
    Xf = X;
    return;
end

% Try toolbox first (better frequency response)
if exist('butter', 'file') == 2 || exist('butter', 'builtin') == 5
    try
        [b, a] = butter(2, fc / (fs / 2), 'low');
        Xf     = zeros(size(X));
        for k = 1:3
            Xf(:, k) = filtfilt(b, a, X(:, k));
        end
        return;
    catch
        % fall through to moving average
    end
end

% Moving-average fallback (symmetric → zero-phase, ~box LP filter)
w  = max(3, round(fs / fc));
Xf = zeros(size(X));
h  = ones(w, 1) / w;
for k = 1:3
    col    = X(:, k);
    % Apply twice in opposite directions for zero-phase
    fwd    = conv(col,  h, 'same');
    rev    = conv(fwd(end:-1:1), h, 'same');
    Xf(:,k) = rev(end:-1:1);
end
end

function m = yn_robust_mean(X, trim_pct)
% Trimmed mean per column (no Statistics Toolbox required).
m = zeros(3, 1);
for k = 1:3
    col = sort(X(:, k));
    n   = numel(col);
    lo  = floor(trim_pct / 100 * n) + 1;
    hi  = ceil((1 - trim_pct / 100) * n);
    lo  = min(lo, hi);
    m(k) = mean(col(lo:hi));
end
end

function s = yn_setdef(s, field, value)
if ~isfield(s, field) || isempty(s.(field))
    s.(field) = value;
end
end
