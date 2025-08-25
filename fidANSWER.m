clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model answer for simulating a free induction decay on an 
%%% electron spin-1/2.
%%%
%%%
%%% Gabriel Moise, Oxford, 2025
%%% https://github.com/moisegabriel/Spin-Dynamics-Without-Tears
%%% αβγδεζηθικλμνξοπρςστυφχψω
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I. USER INPUTS.
% I.1 Spin system parameters: g-factor.
g = 2.00232;

% I.2 Experimental parameters for an FID experiment
B0 = 350; % mT, static magnetic field
mw = 9.808681559288270; % GHz, microwave frequency
B1 = 0.5575; % mT, pulse strength
tp = 16; % ns, pulse length
tau = 300; % ns, free evolution period
dt = 0.25; % ns, time increment (for propagation)

% I.3 Detection operators
detOps = {'Sx', 'Sy', 'Sz'};

% II. PHYSICAL CONSTANTS AND SPIN OPERATORS.
% II.0 Bohr magneton and Planck's constants (SI units)
bmagn = 9.2740e-24; % J*T-1, Bohr magneton
planck = 6.6261e-34; % J*s, Planck's constant
% II.1 Spin operators
E = eye(2); Sz = [0.5 0; 0 -0.5];
Sx = [0 1/2; 1/2 0]; Sy = [0 -1i/2; 1i/2 0];

% III. STATIC HAMILTONIAN IN ROTATING FRAME. 
% Type the formula for the Hamiltonian of an electron in a static 
% magnetic field along the z-axis. The units of H0 should be rad*MHz, 
% so think about how to combine the values of B0, bmagn, and planck.
% Finally, consult the Worksheet for the conversion to the rotating frame.
omega0 = 2*pi*(bmagn*B0*g/planck)*1e-9; % radians*MHz
omegaMW = 2*pi*mw*1e3; % radians*MHz
H0 = (omega0 - omegaMW)*Sz;

% IV. TOTAL HAMILTONIAN DURING PULSE. 
% Using the same logic as above, construct the total Hamiltonian 
% during the microwave pulse.
omega1 = 2*pi*(bmagn*B1*g/planck)*1e-9; % radians*MHz
Hpulse = H0 + omega1*Sx;

% V. INITIAL STATE.
% Consult exercise 7, from Section 2 of the worksheet.
RHO = -Sz;

% VI. DETECTION OPERATORS. 
% In section I.3, the detection operators are defined as strings. 
% Use the template switch statement below to construct a 
% cell of observables which contains matrices instead of strings. 
for iOp = 1:numel(detOps)
    switch detOps{iOp}
        case 'Sz'
            observables{iOp} = Sz;
        case 'Sy'
            observables{iOp} = Sy;
        case 'Sx'
            observables{iOp} = Sx;
    end
end

% VII. PROPAGATION AND DETECTION.
% VII.0 Detect initial signal
for iOp = 1:numel(observables)
    signal{iOp}(1) = trace(RHO*observables{iOp}); 
end

% VII.1 Propagation during the pulse.
% Step 1. Time axis during pulse. The initial time point is dt, the 
% increment is dt, and the final time point is ?
time1 = dt : dt : tp; 
nT1 = length(time1);
% Step 2. Calculate the propagator (U) which takes the density matrix 
% from some initial time t to the next time point (t+dt). The conjugate
% transpose of the propagator is Ud. 
U = expm(-1i*Hpulse*dt*1e-3);
Ud = U';
% Step 3. Loop over each time point. There will be nT time points. 
for iT = 1:nT1
    % Step 4. Propagate! In other words, update the density matrix such
    % that it gets replaced/overwritten by the next density matrix.
    RHO = U*RHO*Ud;
    % Step 5. Detect! Loop over each element of the observables cell and
    % calculate the signal at the current time point. Remember: you
    % already calculated the initial signal!
    for iOp = 1:numel(observables)
        signal{iOp}(1 + iT) = trace(RHO*observables{iOp}); 
    end
end

% VII.2 Propagation during free evolution. Repeat the above procedure for
% the free evolution period.
% Step 1. Time axis during free evolution. The initial time point is ?, the 
% increment is dt, and the final time point is ?+?
time2 = (tp + dt) : dt : (tp + tau);
nT2 = length(time2);
% Step 2. Re-calculate the propagator.
U = expm(-1i*H0*dt*1e-3);
Ud = U';
for iT = 1:nT2
    % Step 3. Propagate!
    RHO = U*RHO*Ud;
    % Step 4: Detect!
    for iOp = 1:numel(observables)
        signal{iOp}(nT1 + 1 + iT) = trace(RHO*observables{iOp});
    end
end

% VIII. PLOT RESULTS.
% Step 1. Assemble the final time axis.
time = [0, time1, time2];
% Step 2. Create a figure and loop over the observables.
figure
hold on; axis tight; box on;
for iOp = 1:numel(observables)
    plot(time,real(signal{iOp}),'LineWidth',1.5,'DisplayName',detOps{iOp})
end
ylabel('Signal Intensity')
xlabel('time / ns')
legend show;