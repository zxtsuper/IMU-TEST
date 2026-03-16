function imu = read_imu_csv(filepath, opts)
% READ_IMU_CSV  Read a 7-column comma-separated IMU data file.
%
%   IMU = READ_IMU_CSV(FILEPATH) reads the file at FILEPATH and returns a
%   struct with fields:
%       t     - Nx1 timestamp vector (seconds)
%       acc   - Nx3 accelerometer matrix [ax, ay, az]  (m/s^2)
%       gyro  - Nx3 gyroscope matrix    [gx, gy, gz]  (rad/s)
%       dt    - median sample interval  (s)
%       fs    - nominal sample rate     (Hz)
%       N     - number of samples
%
%   IMU = READ_IMU_CSV(FILEPATH, OPTS) allows optional configuration via a
%   struct OPTS with fields:
%       gyro_unit  - 'rad_s' (default) or 'deg_s'
%                   Set to 'deg_s' if gyro columns are in deg/s; the
%                   function will convert to rad/s automatically.
%       acc_unit   - 'm_s2' (default) or 'g'
%                   Set to 'g' if acc columns are in gravitational units;
%                   the function will convert to m/s^2 automatically.
%       verbose    - true (default) or false
%
%   File format (no header line):
%       timestamp, ax, ay, az, gx, gy, gz
%   where timestamp is a Unix epoch float, ax/ay/az in m/s^2 and
%   gx/gy/gz in rad/s (default units).
%
%   Example:
%       imu = read_imu_csv('../data/test.txt');
%       plot(imu.t - imu.t(1), imu.acc(:,3));

% ----- defaults --------------------------------------------------------
if nargin < 2 || isempty(opts)
    opts = struct();
end
opts = imu_setdefault(opts, 'gyro_unit', 'rad_s');
opts = imu_setdefault(opts, 'acc_unit',  'm_s2');
opts = imu_setdefault(opts, 'verbose',   true);

% ----- read file -------------------------------------------------------
if ~exist(filepath, 'file')
    error('read_imu_csv: file not found: %s', filepath);
end

raw = dlmread(filepath, ','); %#ok<DLMRD>  % works in R2020

if size(raw, 2) < 7
    error('read_imu_csv: expected >= 7 columns, got %d.', size(raw, 2));
end

t    = raw(:, 1);
acc  = raw(:, 2:4);
gyro = raw(:, 5:7);

% ----- unit conversion ------------------------------------------------
g0 = 9.80665;
if strcmpi(opts.acc_unit, 'g')
    acc = acc * g0;
end
if strcmpi(opts.gyro_unit, 'deg_s')
    gyro = gyro * (pi / 180);
end

% ----- timing ---------------------------------------------------------
dt_vec = diff(t);
dt     = median(dt_vec);
fs     = 1 / dt;

if opts.verbose
    fprintf('[read_imu_csv] Loaded %d samples | fs=%.1f Hz | duration=%.2f s\n', ...
        numel(t), fs, t(end)-t(1));
end

% ----- output struct --------------------------------------------------
imu.t    = t;
imu.acc  = acc;
imu.gyro = gyro;
imu.dt   = dt;
imu.fs   = fs;
imu.N    = numel(t);
end

% ======================================================================
function s = imu_setdefault(s, field, value)
if ~isfield(s, field) || isempty(s.(field))
    s.(field) = value;
end
end
