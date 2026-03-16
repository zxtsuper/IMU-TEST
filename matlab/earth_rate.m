function w_ie = earth_rate(lat_deg, frame)
% EARTH_RATE  Return the Earth rotation rate vector for a given latitude.
%
%   W_IE = EARTH_RATE(LAT_DEG) returns the Earth rotation rate vector in
%   the ENU (East-North-Up) navigation frame as a 3x1 column vector
%   [rad/s]:
%
%       W_IE = [0; Omega*cos(lat); Omega*sin(lat)]
%
%   where Omega = 7.292115e-5 rad/s is the Earth's sidereal rotation rate.
%
%   W_IE = EARTH_RATE(LAT_DEG, FRAME) returns the vector in the requested
%   frame. Supported frames:
%       'ENU' (default) - East-North-Up  [0; Om*cos; Om*sin]
%       'NED'           - North-East-Down [Om*cos; 0; -Om*sin]
%
%   Input:
%       LAT_DEG - geodetic latitude in degrees (positive North)
%                 Default: 31.23 deg N (Shanghai) if not provided.
%       FRAME   - coordinate frame string (default 'ENU')
%
%   Example:
%       w = earth_rate(39.9);   % Beijing
%       fprintf('|w_ie| = %.4e rad/s (should be ~7.29e-5)\n', norm(w));

if nargin < 1 || isempty(lat_deg)
    lat_deg = 31.23;   % default: Shanghai, China
    warning('earth_rate: latitude not specified. Using default %.2f deg N.', lat_deg);
end
if nargin < 2 || isempty(frame)
    frame = 'ENU';
end

Omega = 7.292115e-5;   % Earth sidereal rotation rate [rad/s]
lat   = lat_deg * (pi / 180);

switch upper(frame)
    case 'ENU'
        w_ie = [0; Omega * cos(lat); Omega * sin(lat)];
    case 'NED'
        w_ie = [Omega * cos(lat); 0; -Omega * sin(lat)];
    otherwise
        error('earth_rate: unknown frame ''%s''. Use ''ENU'' or ''NED''.', frame);
end
end
