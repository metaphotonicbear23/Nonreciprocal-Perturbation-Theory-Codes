clear; clc;

%% ================================================================
%  First-order perturbation frequency shift for magneto-optic InAs
%  Fields assumed:
%       E in V/m
%       H in A/m
%       x,y coordinates in meters
%
%  Output:
%       delta_omega in rad/s
%       delta_f in THz
%% ================================================================


%% ---------------- Physical constants ----------------------------
eps0 = 8.854187817e-12;      % vacuum permittivity [F/m]
mu0  = 4*pi*1e-7;            % vacuum permeability [H/m]

%% ---------------- Material parameters ---------------------------
eps_inf = 12.3;              % high-frequency permittivity of InAs

% Plasma frequency [Hz]
omega_p = 2.7422e14/(2*pi);

% Cyclotron frequency
%Actual 0.05 T [Hz]
omega_c = -0.5*5.3347e11/(2*pi);

%% ---------------- Resonance frequency ---------------------------
%Lossy
Gamma = 1.55e11;

%Lossless
%Gamma = 0;

%% ---------------- Resonance frequency ---------------------------
% Zero-field resonance angular frequency [Hz]
%structure 1`
%omega = 1.224789744427E+13;

%structure 1 InGaAs on InAs

%67deg
%omega =  1.2300691695388E+13;

%72deg
%omega = 1.227596366315E+13;

%73deg
%omega = 12262453288613;

%74deg
%omega = 1.2257439610761E+13;pert

%75deg
%omega = 1.2252430031061E+13;

%76deg
%omega = 1.2247424544489E+13;

%77deg
%omega = 1.2243423098914E+13;

%78deg
%omega = 1.2239424267167E+13;

%79deg
omega = 1.2235428046690E+13;

%% ---------------- Load field profiles ---------------------------
% Expected columns:
%   column 1: x-coordinate [m]
%   column 2: y-coordinate or z-coordinate [m]
%   column 3: field value

data_Ex_real = readmatrix('real_Ex_79deg_24pt502um_0T.txt');
data_Ex_imag = readmatrix('imag_Ex_79deg_24pt502um_0T.txt');

data_Ey_real = readmatrix('real_Ey_79deg_24pt502um_0T.txt');
data_Ey_imag = readmatrix('imag_Ey_79deg_24pt502um_0T.txt');

data_Hz_real = readmatrix('real_Hz_79deg_24pt502um_0T.txt');
data_Hz_imag = readmatrix('imag_Hz_79deg_24pt502um_0T.txt');

%% ---------------- Extract coordinates ---------------------------
x = data_Ex_real(:,1);       % [m]
y = data_Ex_real(:,2);       % [m]

%% ---------------- Construct complex fields ----------------------
Ex = data_Ex_real(:,3) + 1i*data_Ex_imag(:,3);
Ey = data_Ey_real(:,3) + 1i*data_Ey_imag(:,3);
Hz = data_Hz_real(:,3) + 1i*data_Hz_imag(:,3);

%% ---------------- Basic consistency checks ----------------------
N = length(Ex);

if length(Ey) ~= N || length(Hz) ~= N
    error('Ex, Ey, and Hz arrays do not have the same length.');
end

if any(abs(data_Ey_real(:,1) - x) > 1e-15) || any(abs(data_Hz_real(:,1) - x) > 1e-15)
    warning('x-coordinates may not match between field files.');
end

if any(abs(data_Ey_real(:,2) - y) > 1e-15) || any(abs(data_Hz_real(:,2) - y) > 1e-15)
    warning('y-coordinates may not match between field files.');
end

%% ---------------- Build integration weights ---------------------
% Important:
% Use coordinates in meters for physical SI integration.
% Do NOT multiply by 1e9 here unless you intentionally want nm^2 weights.

pts = [x, y];

tri = delaunayTriangulation(pts);

weights = compute_triangle_area_weights(tri, pts);   % units: m^2

%% ================================================================
%  Denominator: total modal energy per unit length
%
%  For lossless Drude:
%
%      epsilon(omega) = eps_inf * (1 - omega_p^2/omega^2)
%
%  dispersive energy contribution:
%
%      1/2 eps0 eps_inf (omega_p^2/omega^2) |E|^2
%
%  so denominator density:
%
%      W = 1/2 eps0 eps_inf |E|^2
%        + 1/2 mu0 |H|^2
%        + 1/2 eps0 eps_inf (omega_p^2/omega^2) |E|^2
%
%% ================================================================

E2 = abs(Ex).^2 + abs(Ey).^2;
H2 = abs(Hz).^2;

E_energy = 0.5 * eps0 * eps_inf .* E2;

H_energy = 0.5 * mu0 .* H2;

Drude_energy = 0.5 * eps0 * eps_inf * (omega_p^2/omega^2) .* E2;

W_density = E_energy + H_energy + Drude_energy;

W_integral = sum(W_density .* weights);

%% ================================================================
%  Numerator:
%
%      -eps0 * (omega_c omega_p^2 eps_inf / omega^2)
%       * Im(Ex^* Ey - Ey^* Ex)
%
%  Note:
%
%      Im(Ex^* Ey - Ey^* Ex) = 2 Im(Ex^* Ey)
%
%% ================================================================

spin_E = imag(conj(Ex).*Ey - conj(Ey).*Ex);

lossy = (omega_p^4 * Gamma^2 * eps_inf^2 + omega^2 * omega_p^4 * eps_inf^2)/(omega^2+Gamma^2)^2;

num_density = -eps0 * (omega_c /(omega_p^2 * eps_inf) ) *lossy.* spin_E;

num_integral = sum(num_density .* weights);


%% ---------------- Frequency shift -------------------------------
delta_omega = num_integral / W_integral;      % [rad/s]
delta_f_Hz  = delta_omega / (2*pi);           % [Hz]
delta_f_THz = delta_f_Hz / 1e12;              % [THz]

%% ---------------- Print results --------------------------------
fprintf('\n================ Perturbation Result ================\n');
fprintf('Numerator integral:     %.6e\n', num_integral);
fprintf('Denominator integral:   %.6e\n', W_integral);
fprintf('Delta omega:            %.6e rad/s\n', delta_omega);
fprintf('Delta f:                %.6e Hz\n', delta_f_Hz);
fprintf('Delta f:                %.6e THz\n', delta_f_THz);
fprintf('=====================================================\n\n');


%% ---------------- Print energy diagnostics ----------------------
E_int      = sum(E_energy .* weights);
H_int      = sum(H_energy .* weights);
Drude_int  = sum(Drude_energy .* weights);

fprintf('Energy contributions:\n');
fprintf('Electric background:    %.6e\n', E_int);
fprintf('Magnetic:               %.6e\n', H_int);
fprintf('Drude kinetic:          %.6e\n', Drude_int);
fprintf('Total denominator:      %.6e\n\n', W_integral);

fprintf('Energy fractions:\n');
fprintf('Electric fraction:      %.4f\n', E_int/W_integral);
fprintf('Magnetic fraction:      %.4f\n', H_int/W_integral);
fprintf('Drude fraction:         %.4f\n\n', Drude_int/W_integral);


%% ================================================================
%  Local function: triangle-area integration weights
%% ================================================================
function weights = compute_triangle_area_weights(tri, pts)

    N = size(pts, 1);
    weights = zeros(N, 1);

    T = tri.ConnectivityList;
    P = tri.Points;

    for ii = 1:size(T, 1)

        vertex_indices = T(ii, :);
        vertices = P(vertex_indices, :);

        area = polyarea(vertices(:,1), vertices(:,2));

        % Distribute triangle area equally to its three vertices
        weights(vertex_indices) = weights(vertex_indices) + area/3;

    end

end