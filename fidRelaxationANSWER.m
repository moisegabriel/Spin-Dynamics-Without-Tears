clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulating a free induction decay on an electron spin-1/2.
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
T1 = 10; % μs
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
E_ = reshape(E,4,[]); Sz_ = reshape(Sz,4,[]);
Sx_ = reshape(Sx,4,[]); Sy_ = reshape(Sy,4,[]);

% III. STATIC HAMILTONIAN IN ROTATING FRAME
omega0 = 2*pi*(bmagn*B0*g/planck)*1e-9; % radians*MHz
omegaMW = 2*pi*mw*1e3; % radians*MHz
% III.1 Hamiltonian operator in Hilbert space
H0 = (omega0 - omegaMW)*Sz;
% III.2 Hamiltonian commutation superoperator
H0super = kron(E,H0) - kron(transpose(H0),E);

% IV. PULSE HAMILTONIAN IN ROTATING FRAME
omega1 = 2*pi*(bmagn*B1*g/planck)*1e-9; % radians*MHz
% IV.1 Hamiltonian operator in Hilbert space
Hpulse = H0 + omega1*Sx;
% IV.2 Hamiltonian commutation superoperator
Hpulsesuper = kron(E,Hpulse) - kron(transpose(Hpulse),E);

% V. EQUILIBRIUM STATE
EQ = -Sz_;

% VI. RELAXATION SUPEROPEARTOR
Gamma = [0 0 0 -1/T1;
         0 1/T2 0 0;
         0 0 1/T2 0;
         -1/T1 0 0 0];

% VII. OBSERVABLE VECTORS IN LIOUVILLE SPACE
for iOp = 1:numel(detOps)
    switch detOps{iOp}
        case 'x'
            observables{iOp} = Sx_;
        case 'y'
            observables{iOp} = Sy_;
        case 'z'
            observables{iOp} = Sz_;
    end
end

% VIII. PROPAGATION AND DETECTION
% VIII.0 Detect initial signal
RHO = EQ;
for iOp = 1:numel(observables)
    signal{iOp}(1) = RHO'*observables{iOp}; 
end
% VIII.1 Propagation during the pulse
time1 = dt:dt:tp;
nT = length(time1);
L = -1i*Hpulsesuper - Gamma;
SS = -L\(Gamma*EQ);
U = expm(L*dt*1e-3);
for iT = 1:nT
    % Propagate
    RHO = SS + U*(RHO-SS);
    % Detect
    for iOp = 1:numel(observables)
        signal{iOp}(1 + iT) = RHO'*observables{iOp}; 
    end
end
counter = 1 + nT;

% VIII.2 Propagation during free evolution
time2 = (tp+dt):dt:(tp+tau);
nT = length(time2);
L = -1i*H0super - Gamma;
SS = -L\(Gamma*EQ);
U = expm(L*dt*1e-3);
for iT = 1:nT
    % Propagate
    RHO = SS + U*(RHO-SS);
    % Detect
    for iOp = 1:numel(observables)
        signal{iOp}(counter + iT) = RHO'*observables{iOp}; 
    end
end

% IX. PLOT RESULTS
time = [0, time1, time2];
figure
hold on; axis tight; box on;
for iOp = 1:numel(observables)
    plot(time,real(signal{iOp}),'LineWidth',1.5,'DisplayName',detOps{iOp})
end
xline(tp,'r--','HandleVisibility','off')
xlabel('time / ns')
legend show;


            
