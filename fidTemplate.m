clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Template for simulating a free induction decay on an electron spin-1/2.
%%%
%%%
%%%
%%% Gabriel Moise, Oxford, 2025
%%% αβγδεζηθικλμνξοπρςστυφχψω
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I. USER INPUTS
% I.1 Spin system parameters: g-factor.
g = 2.00232;

% I.2 Experimental parameters. Consider the pulse sequence for the FID
% experiment. What are the experimental parameters? Define them here as
% variables. Make sure to specify the units and choose sensible 
% variable names.
B0 = 350; % mT, static field

% I.3 Detection operators
detOps = {'Sx', 'Sy', 'Sz'};

% II. PHYSICAL CONSTANTS AND SPIN OPERATORS
% II.0 Bohr magneton and Planck's constants (SI units)
bmagn = 9.2740e-24; % J*T-1
planck = 6.6261e-34; % J*s
% II.1 Spin operators
E = eye(2); Sz = [0.5 0; 0 -0.5];
Sx = ?; Sy = ?;

% III. STATIC HAMILTONIAN IN ROTATING FRAME
H0 = ?;

% IV. PULSE HAMILTONIAN IN ROTATING FRAME
Hpulse = H0 + ?;

% V. INITIAL STATE
RHO = ?;

% VI. DETECTION OPERATORS. Add the remaining observables.
for iOp = 1:numel(detOps)
    switch detOps{iOp}
        case 'Sz'
            observables{iOp} = Sz;
    end
end

% VII. PROPAGATION AND DETECTION
% VII.0 Detect initial signal
for iOp = 1:numel(observables)
    signal{iOp}(1) = ?; 
end
% VII.1 Propagation during the pulse

U = ?;
Ud = ?;
% Loop over each time point. There will be nT time points. 
for iT = 1:nT
    % Propagate
    RHO = ?;
    % Detect
    for iOp = 1:numel(observables)
        signal{iOp}(?) = ?; 
    end
end

% VII.2 Propagation during free evolution

U = ?;
Ud = ?;
for iT = 1:nT
    % Propagate
    RHO = ?;
    % Detect
    for iOp = 1:numel(observables)
        signal{iOp}(?) = ?; 
    end
end

% VIII. PLOT RESULTS
time = ?;
figure
hold on; axis tight; box on;
for iOp = 1:numel(observables)
    plot(time,real(signal{iOp}),'LineWidth',1.5,'DisplayName',detOps{iOp})
end
ylabel('Signal Intensity / ?')
xlabel('time / ns')
legend show;