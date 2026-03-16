# IMU-TEST — MEMS-IMU Quick North-Finding (Gyrocompass) Demo

MATLAB R2020 example project for estimating compass heading (yaw / north-finding) using **only a MEMS IMU** (accelerometer + gyroscope). No GNSS, magnetometer, or other external sensors are used.

---

## Contents

```
IMU-TEST/
├── data/
│   └── test.txt          – sample IMU data (30 s, 200 Hz, synthetic)
└── matlab/
    ├── north_finding_demo.m           – main runnable script  ← start here
    ├── read_imu_csv.m                 – 7-column CSV data reader
    ├── detect_static_segments.m      – automatic static window detection
    ├── earth_rate.m                   – Earth rotation rate vector (ENU/NED)
    └── estimate_yaw_northfinding.m   – EKF gyrocompass core solver
```

---

## Quick Start

1. **Open MATLAB R2020** (or later).
2. Set your geodetic latitude in `north_finding_demo.m`:
   ```matlab
   LATITUDE_DEG = 31.23;   % degrees North; e.g. 39.9 = Beijing
   ```
3. Point `DATA_FILE` at your IMU file (default `../data/test.txt`).
4. Run:
   ```matlab
   cd matlab
   north_finding_demo
   ```

The script outputs:
- Estimated **heading** (compass bearing of IMU x-axis, degrees from North, clockwise-positive).
- 1-sigma uncertainty estimate.
- Roll and pitch (from accelerometer gravity alignment).
- Four diagnostic figures.

---

## Input Data Format

Plain-text, comma-separated, **no header line**, 7 columns per row:

```
timestamp, ax, ay, az, gx, gy, gz
```

| Column | Meaning | Default unit |
|--------|---------|--------------|
| 1 | Unix epoch timestamp (float seconds) | s |
| 2 | Accelerometer x | m/s² |
| 3 | Accelerometer y | m/s² |
| 4 | Accelerometer z | m/s² |
| 5 | Gyroscope x | rad/s |
| 6 | Gyroscope y | rad/s |
| 7 | Gyroscope z | rad/s |

Example rows:
```
1716340748.057519, 0.008340, 0.003687, -9.807658, 0.000203, 0.000203, 0.000030
1716340748.062519, -0.000755, 0.004580, -9.800848, 0.000121, 0.000168, 0.000011
```

If your gyroscope data is in **deg/s**, set `GYRO_UNIT = 'deg_s'` in the demo script. Similarly use `ACC_UNIT = 'g'` for gravitational-unit accelerometers.

---

## Algorithm Overview

### 1 — Static Window Detection
The script automatically finds the quietest contiguous segment of data where:
- The accelerometer magnitude is close to *g* (device is not moving).
- The gyroscope standard deviation is minimal (no rotation noise).

This avoids using noisy or dynamic segments for alignment.

### 2 — Coarse Alignment (roll & pitch)
With the sensor stationary, the mean accelerometer reading is a projection of the gravity vector onto the body frame. From this, roll and pitch are computed analytically:
```
pitch = atan2(ux, sqrt(uy² + uz²))
roll  = atan2(-uy, uz)
where u = normalize(mean_acc)
```

### 3 — EKF Gyrocompass (heading / yaw estimation)

**State** (4 × 1):
```
x = [ψ; b_gx; b_gy; b_gz]
```
where `ψ` is the yaw/heading angle and `b_g` is the constant gyroscope bias vector.

**Process model** (static, random-walk):
```
ψ_{k+1}  = ψ_k  + w_ψ
b_{k+1}  = b_k  + w_b
```

**Measurement model**:
```
z_k = gyro_body
h(x) = Cbn(ψ, roll, pitch) · ω_ie_ENU + b_g
```
where `Cbn = Rx(φ) · Ry(θ) · Rz(ψ)` rotates the known Earth rotation rate vector from the ENU navigation frame into the body frame. The Jacobian is computed numerically.

The filter simultaneously estimates the heading **and** cancels the gyroscope bias, making it robust to the 5°/h drift level typical of MEMS sensors.

**Noise Suppression**:
- Zero-phase low-pass filter on acc/gyro before the EKF (default 1.5 Hz cut-off, uses `butter`+`filtfilt` if Signal Processing Toolbox is present, otherwise a symmetric moving average).
- Trimmed mean for roll/pitch update inside the filter loop.
- Automatic outlier rejection via static-segment selection.

### 4 — Heading Output Convention
```
heading_deg = atan2d(cos(ψ), -sin(ψ))  [0–360, North = 0, clockwise+]
```
This is the compass bearing of the **IMU x-axis** (ax column). To obtain the bearing of a different body axis, apply the appropriate offset.

---

## Coordinate Convention

| Frame | Definition |
|-------|-----------|
| Navigation | **ENU**: x = East, y = North, z = Up |
| Body | Generic; body z typically points **down** (az ≈ −g at rest) |
| `Cbn` | `Rx(φ) · Ry(θ) · Rz(ψ)` — rotates ENU vectors into body frame |

At rest: `acc_body ≈ Cbn · [0; 0; g₀]`  
Gyro measurement: `gyro_body ≈ Cbn · ω_ie_ENU + bias`

---

## Configuring the Latitude

Latitude must be set to your actual location for accurate north-finding. The Earth rotation rate horizontal component (`Ω·cos(lat)`) is the primary observable signal; a wrong latitude will introduce a systematic heading error.

```matlab
LATITUDE_DEG = 31.23;   % Shanghai
LATITUDE_DEG = 39.91;   % Beijing
LATITUDE_DEG = 22.52;   % Guangzhou
LATITUDE_DEG = 51.51;   % London
```

---

## Tuning Parameters for Your Sensor

| Parameter | Default | Effect |
|-----------|---------|--------|
| `MAX_EKF_S` | 60 s | Longer → better accuracy, more temperature drift risk. Cap at 120–180 s. |
| `MIN_STATIC_S` | 10 s | Minimum data required. Increase if noisy environment. |
| `LP_CUTOFF_HZ` | 1.5 Hz | Lower → smoother but more lag. 0.5–3 Hz is typical for static. |
| `SIGMA_GYRO` | 2e-5 rad/s | **Post-lowpass** measurement noise std. Set ≈ std(filtered gyro static). |
| `MAX_BIAS` | 3e-4 rad/s | Expected max gyro bias (1-sigma prior). For 5°/h: ~2.4e-4 rad/s. |
| `SIGMA_BIAS_RW` | 2e-6 rad/s/√s | Larger allows faster bias tracking; too large → noisy estimate. |

**How to measure your `SIGMA_GYRO`:**
```matlab
imu = read_imu_csv('../data/test.txt');
% after low-pass filter — see estimate_yaw_northfinding.m for yn_lowpass
fprintf('gyro x std = %.3e rad/s\n', std(imu.gyro(:,1)));
```
Use the per-axis standard deviation of a quiet static segment as `SIGMA_GYRO`.

---

## Expected Accuracy (5°/h MEMS Gyro)

Earth rotation rate magnitude ≈ 7.3×10⁻⁵ rad/s ≈ 15°/h.  
A 5°/h gyro bias is **one-third of the earth-rate signal**; north-finding is challenging but feasible.

| Static duration | Typical heading accuracy (1σ) |
|-----------------|-------------------------------|
| 30 s | 15 – 30° |
| 60 s | 5 – 15° |
| 120 s | 3 – 8° |
| 180 s | 2 – 6° |

> **Practical tips:**
> - Place the sensor on a vibration-free surface (not a running machine).
> - Avoid air drafts and cable tension that can cause micro-vibrations.
> - If the gyro white-noise floor is lower than the 5°/h bias (common in higher-grade MEMS), accuracy improves significantly.
> - For latitudes below ~15°N/S the horizontal Earth-rate component is very small and accuracy degrades rapidly.

---

## Dependencies

| Requirement | Notes |
|-------------|-------|
| MATLAB R2020 | Core functionality; no newer syntax used |
| Signal Processing Toolbox | **Optional**. Used for `butter`+`filtfilt`; falls back to moving-average if absent |
| Statistics Toolbox | **Not required**. Trimmed mean is implemented manually |
| Optimization Toolbox | **Not required**. EKF uses only matrix arithmetic |

---

## File Format Reference

See `matlab/read_imu_csv.m` for the full parser. To read data in code:
```matlab
imu = read_imu_csv('../data/test.txt');
% imu.t     Nx1 timestamps (s)
% imu.acc   Nx3 accelerometer (m/s^2)
% imu.gyro  Nx3 gyroscope (rad/s)
% imu.fs    sample rate (Hz)
```