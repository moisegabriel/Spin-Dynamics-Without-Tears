clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulating pulse EPR sequences on an electron spin-1/2
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
B1 = 0.55; % mT, pulse strength
tp = 16; % ns, pulse length
tau = 300; % ns, free evolution period
dt = 0.25; % ns, time increment
sequence = {{'x',tp,B1},{'free',tau},{'x',2*tp,B1},{'free',tau+2*tp}};
detOps = {'Sx', 'Sy', 'Sz'};

% II. PHYSICAL CONSTANTS AND SPIN OPERATORS
bmagn = 9.2740e-24; % J*T-1
planck = 6.6261e-34; % J*s
E = eye(2); Sz = [0.5 0; 0 -0.5];
Sx = [0 0.5; 0.5 0]; Sy = [0 -0.5i; 0.5i 0];

% III. STATIC HAMILTONIAN IN ROTATING FRAME
omega0 = 2*pi*(bmagn*B0*g/planck)*1e-9; % radians*MHz
omegaMW = 2*pi*mw*1e3; % radians*MHz
H0 = omega0*Sz - omegaMW*Sz;

% IV. CONSTRUCT PROPAGATORS AND TIME AXIS
time = [0];
for iS = 1:numel(sequence)
    currentTime = (time(end)+dt):dt:(time(end) + sequence{iS}{2});
    time = [time, currentTime];
    nT{iS} = length(currentTime);
    switch sequence{iS}{1}
        case 'free'
            Propagator{iS} = expm(-1i*H0*dt*1e-3);
        case 'x'
            omega1 = 2*pi*(bmagn*sequence{iS}{3}*g/planck)*1e-9; % radians*MHz
            Propagator{iS} = expm(-1i*(H0 + omega1*Sx)*dt*1e-3);
    end
end

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
% VII.1 Loop over sequence
counter = 1;
for iS = 1:length(sequence)
    % Retrieve propagator for current event
    U = Propagator{iS};
    Ud = U';
    % Loop over the time points of the current event
    for iT = 1:nT{iS}
        counter = counter + 1;
        % Propagate
        RHO = U*RHO*Ud; 
        % Detect
        for iOp = 1:numel(observables)
            signal{iOp}(counter) = trace(RHO*observables{iOp}); 
        end
        
    end
end

% VIII. PLOT RESULTS
figure
hold on; axis tight; box on;
for iOp = 1:numel(observables)
    plot(time,real(signal{iOp}),'LineWidth',1.5,'DisplayName',detOps{iOp})
end
xlabel('time / ns')
xline(tp,'r--','HandleVisibility','off')
xline(tau+tp,'r--','HandleVisibility','off')
xline(tau+3*tp,'r--','HandleVisibility','off')
legend show;