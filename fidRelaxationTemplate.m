clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Template for simulating a free induction decay on an electron spin-1/2.
%%% Including relaxation.
%%%
%%%
%%%
%%% Gabriel Moise, Oxford, 2025
%%% αβγδεζηθικλμνξοπρςστυφχψω
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I. USER INPUTS
% I.1 Spin system parameters
g = 2.00232;
T1 = 1e10; % μs
T2 = 0.5; % μs

% I.2 Experimental parameters
B0 = 350; % mT, static field
mw = 9.808681559288270; % GHz, microwave frequency
B1 = 0.5575; % mT, pulse strength
tp = 16; % ns, pulse length
tau = 300; % ns, free evolution period
dt = 0.25; % ns, time increment
detOps = {'x', 'y', 'z'};

% II. PHYSICAL CONSTANTS AND SPIN OPERATORS
bmagn = 9.2740e-24; % J*T-1
planck = 6.6261e-34; % J*s
% II.1 Operators as 2-by-2 matrices acting on the Hilbert space
E = eye(2); Sz = [0.5 0; 0 -0.5];
Sx = [0 0.5; 0.5 0]; Sy = [0 -0.5i; 0.5i 0];
% II.2 Promotion to vectors in Liouville space
E_ = reshape(E,?,?); Sz_ = ?;
Sx_ = ?; Sy_ = ?;

% III. STATIC HAMILTONIAN IN ROTATING FRAME
omega0 = 2*pi*(bmagn*B0*g/planck)*1e-9; % radians*MHz
omegaMW = 2*pi*mw*1e3; % radians*MHz
% III.1 Hamiltonian operator in Hilbert space
H0 = (omega0 - omegaMW)*Sz;
% III.2 Hamiltonian commutation superoperator
H0super = ?;

% IV. PULSE HAMILTONIAN IN ROTATING FRAME
omega1 = 2*pi*(bmagn*B1*g/planck)*1e-9; % radians*MHz
% IV.1 Hamiltonian operator in Hilbert space
Hpulse = H0 + omega1*Sx;
% IV.2 Hamiltonian commutation superoperator
Hpulsesuper = ?;

% V. EQUILIBRIUM STATE
EQ = ?;

% VI. RELAXATION SUPEROPEARTOR
Gamma = ?;

% VII. OBSERVABLE VECTORS IN LIOUVILLE SPACE
for iOp = 1:numel(detOps)
    switch detOps{iOp}
        ?
    end
end

% VIII. PROPAGATION AND DETECTION
% VIII.0 Detect initial signal
RHO = ?;
for iOp = 1:numel(observables)
    signal{iOp}(1) = ?; 
end
% VIII.1 Propagation during the pulse
?
for iT = 1:nT
    % Propagate
    RHO = ?;
    % Detect
    for iOp = 1:numel(observables)
        signal{iOp}(?) = ?; 
    end
end

% VIII.2 Propagation during free evolution
?
for iT = 1:nT
    % Propagate
    RHO = ?;
    % Detect
    for iOp = 1:numel(observables)
        signal{iOp}(?) = ?; 
    end
end

% IX. PLOT RESULTS
time = ?;
figure
hold on; axis tight; box on;
for iOp = 1:numel(observables)
    plot(time,real(signal{iOp}),'LineWidth',1.5,'DisplayName',detOps{iOp})
end
xline(tp,'r--','HandleVisibility','off')
xlabel('time / ns')
legend show;

            