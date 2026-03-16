function [i0, i1, score] = detect_static_segments(t, acc, gyro, opts)
% DETECT_STATIC_SEGMENTS  Find the highest-quality static window in IMU data.
%
%   [I0, I1, SCORE] = DETECT_STATIC_SEGMENTS(T, ACC, GYRO) searches the
%   data for the contiguous window where:
%     - the accelerometer magnitude is closest to g  (device is level/static)
%     - the gyroscope magnitude is near zero         (no rotation)
%     - both have minimal variance                   (quiet environment)
%
%   Inputs:
%       T     - Nx1 timestamp vector (s)
%       ACC   - Nx3 accelerometer matrix (m/s^2)
%       GYRO  - Nx3 gyroscope matrix (rad/s)
%       OPTS  - (optional) configuration struct:
%           min_static_s   - minimum window length in seconds (default 10)
%           max_static_s   - maximum window length in seconds (default 60)
%           gyro_std_th    - gyro magnitude std threshold (default 0.015 rad/s)
%           acc_std_th     - acc magnitude std threshold  (default 0.15 m/s^2)
%           acc_mean_tol   - tolerance |accMag - g| for static (default 0.5 m/s^2)
%           search_step_s  - step size for window search in s (default 0.5 s)
%           verbose        - print summary (default true)
%
%   Outputs:
%       I0    - start index of best static window
%       I1    - end   index of best static window
%       SCORE - quality score (lower is better); Inf if no window found
%
%   If no window passes all thresholds, the function falls back to the
%   first MIN_STATIC_S seconds of data (with a warning).

% ----- defaults -------------------------------------------------------
if nargin < 4 || isempty(opts)
    opts = struct();
end
opts = dss_setdefault(opts, 'min_static_s',  10);
opts = dss_setdefault(opts, 'max_static_s',  60);
opts = dss_setdefault(opts, 'gyro_std_th',   0.015);
opts = dss_setdefault(opts, 'acc_std_th',    0.15);
opts = dss_setdefault(opts, 'acc_mean_tol',  0.5);
opts = dss_setdefault(opts, 'search_step_s', 0.5);
opts = dss_setdefault(opts, 'verbose',       true);

% ----- setup ----------------------------------------------------------
t    = t(:);
N    = numel(t);
dt   = median(diff(t));
fs   = 1 / dt;
g0   = 9.80665;

minL = max(round(opts.min_static_s * fs), 10);
maxL = min(round(opts.max_static_s * fs), N);
step = max(1, round(opts.search_step_s * fs));

acc_mag  = sqrt(sum(acc .^ 2, 2));
gyro_mag = sqrt(sum(gyro .^ 2, 2));

best_score = Inf;
i0 = 1;
i1 = min(minL, N);

% ----- sliding window search -----------------------------------------
for L = minL : step : maxL
    for s = 1 : step : (N - L + 1)
        e = s + L - 1;
        
        a_seg = acc_mag(s:e);
        g_seg = gyro_mag(s:e);
        
        accStd  = std(a_seg);
        gyroStd = std(g_seg);
        accMean = mean(a_seg);
        
        if gyroStd  < opts.gyro_std_th && ...
           accStd   < opts.acc_std_th  && ...
           abs(accMean - g0) < opts.acc_mean_tol
            score = gyroStd + 0.5 * accStd + 0.1 * abs(accMean - g0);
            if score < best_score
                best_score = score;
                i0 = s;
                i1 = e;
            end
        end
    end
end

score = best_score;

if ~isfinite(best_score)
    i1 = min(minL, N);
    i0 = 1;
    score = Inf;
    if opts.verbose
        warning(['detect_static_segments: no window passed thresholds. ' ...
                 'Using first %.1f s. Consider relaxing thresholds or ' ...
                 'ensuring data is from a stationary device.'], ...
                 (i1 - i0) * dt);
    end
else
    if opts.verbose
        fprintf('[detect_static] Best static window: samples %d–%d (%.1f – %.1f s) | %.1f s | score=%.4g\n', ...
            i0, i1, (i0-1)*dt, (i1-1)*dt, (i1-i0)*dt, best_score);
    end
end
end

% ======================================================================
function s = dss_setdefault(s, field, value)
if ~isfield(s, field) || isempty(s.(field))
    s.(field) = value;
end
end
