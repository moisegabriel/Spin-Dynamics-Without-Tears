clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulating a free induction decay on an electron spin-1/2.
%%%
%%%
%%%
%%% Gabriel Moise, Oxford, 2025
%%% αβγδεζηθικλμνξοπρςστυφχψω
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I. USER INPUTS
% I.1 Spin system parameters
g = 2.00232;

% I.2 Experimental parameters
B0 = 350; % mT, static field
mw = 9.808681559288270; % GHz, microwave frequency
B1 = 0.5575; % mT, pulse strength
tp = 16; % ns, pulse length
tau = 300; % ns, free evolution period
dt = 0.25; % ns, time increment
detOps = {'Sx', 'Sy', 'Sz'};

% II. PHYSICAL CONSTANTS AND SPIN OPERATORS
bmagn = 9.2740e-24; % J*T-1
planck = 6.6261e-34; % J*s
E = eye(2); Sz = [0.5 0; 0 -0.5];
Sx = [0 0.5; 0.5 0]; Sy = [0 -0.5i; 0.5i 0];

% III. STATIC HAMILTONIAN IN ROTATING FRAME
omega0 = 2*pi*(bmagn*B0*g/planck)*1e-9; % radians*MHz
omegaMW = 2*pi*mw*1e3; % radians*MHz
H0 = (omega0 - omegaMW)*Sz;

% IV. PULSE HAMILTONIAN IN ROTATING FRAME
omega1 = 2*pi*(bmagn*B1*g/planck)*1e-9; % radians*MHz
Hpulse = H0 + omega1*Sx;

% V. INITIAL STATE
RHO = -Sz;

% VI. DETECTION OPERATORS
for iOp = 1:numel(detOps)
    switch detOps{iOp}
        case 'Sx'
            observables{iOp} = Sx;
        case 'Sy'
            observables{iOp} = Sy;
        case 'Sz'
            observables{iOp} = Sz;
    end
end

% VII. PROPAGATION AND DETECTION
% VII.0 Detect initial signal
for iOp = 1:numel(observables)
    signal{iOp}(1) = trace(RHO*observables{iOp}); 
end
% VII.1 Propagation during the pulse
time1 = dt:dt:tp;
nT = length(time1);
U = expm(-1i*Hpulse*dt*1e-3);
Ud = U';
for iT = 1:nT
    % Propagate
    RHO = U*RHO*Ud;
    % Detect
    for iOp = 1:numel(observables)
        signal{iOp}(1 + iT) = trace(RHO*observables{iOp}); 
    end
end
counter = 1 + nT;

% VII.2 Propagation during free evolution
time2 = (tp+dt):dt:(tp+tau);
nT = length(time2);
U = expm(-1i*H0*dt*1e-3);
Ud = U';
for iT = 1:nT
    % Propagate
    RHO = U*RHO*Ud;
    % Detect
    for iOp = 1:numel(observables)
        signal{iOp}(counter + iT) = trace(RHO*observables{iOp}); 
    end
end

% VIII. PLOT RESULTS
time = [0, time1, time2];
figure
hold on; axis tight; box on;
for iOp = 1:numel(observables)
    plot(time,real(signal{iOp}),'LineWidth',1.5,'DisplayName',detOps{iOp})
end
xline(tp,'r--','HandleVisibility','off')
xlabel('time / ns')
legend show;

            