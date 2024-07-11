% Title: Code to get the fluid flow from a set of Regularized Stokeslets.
% Author: Stephen Williams.

%close all
%clear all %#ok<CLALL>

%% Add the function files need to run
addpath('functions/')
addpath('classes/')

%% Set parameters
parameters % Set the parameters

%% Set channel geometry
stks = getStokesletPositions(rho,system,U0);

%% Solve for the forces
[iS] = getForces(stks,eps_reg);

%% Simulate the particle motion

nparticles = 100; % Number of particle trajectories to simulate
U0 = 100/45; % Pouiselle parameters
delta = 20; % Boundary force potential scaling
tmin = 0; tmax = 15; ntsteps = 100; % Time paramaters
ic = [linspace(-9,9,nparticles);18*ones(1,nparticles)];
omega = 6;
%options = odeset('RelTol', 1e-4);
[t,y] = ode45(@(t,y) odefun(t,y,stks,iS,eps_reg,omega,U0,delta,system),linspace(tmin,tmax,ntsteps),ic);

x_pos = y(:,1:2:end); y_pos = y(:,2:2:end); % Extract the positions into an array for plotting

scatter(stks(:,1),stks(:,2),0.1,'r'); hold on
for i = 1:nparticles
    plot(x_pos(:,i),y_pos(:,i),'k')
end

axis equal

%% Functions needed for ODE45 

function dydt = odefun(t1,y1,stks1,iS1,eps_reg1,omega,U0,delta1,system1)

R = 1; % Radius of appendages
r = 0.2; % Radius of particles

a = find(stks1(:,3) == 8); % Find the Pousielle boundary sections
Ut = U0*(1+cos(t1*omega))/2; % Set the maximum flow rate in
stks1(a,4:5) = poisuelleFlow(length(a),Ut); % Re-set the Pousielle boundary sections
a = find(stks1(:,3) == 9); % Find the Pousielle boundary sections
Ut = U0*(1+cos(t1*omega))/2; % Set the maximum flow rate in
stks1(a,4:5) = poisuelleFlow(length(a),Ut); % Re-set the Pousielle boundary sections

Uflow = calculateFlowVector(stks1,iS1,y1,eps_reg1); % Get the flow from fluid interactions

% Get the appendage centers
AP1 = [system1.appendage_parameters(3) + (1+system1.appendage_parameters(1)/2)*cos(system1.appendage_parameters(2)), ...
       system1.appendage_parameters(4) + (1+system1.appendage_parameters(1)/2)*sin(system1.appendage_parameters(2))];
AP2 = [system1.appendage_parameters(3) - (1+system1.appendage_parameters(1)/2)*cos(system1.appendage_parameters(2)), ... 
       system1.appendage_parameters(4) - (1+system1.appendage_parameters(1)/2)*sin(system1.appendage_parameters(2))];
AP3 = [-1,1].*AP1; AP4 = [-1,1].*AP2;
app_pos = [AP1;AP2;AP3;AP4]; % Collect up all of the the positions
x = [y1(1:2:end)';y1(2:2:end)'];
[~,I] = min([vecnorm(x(1:2,:)-AP1',2) ; ...
             vecnorm(x(1:2,:)-AP2',2) ; ...
             vecnorm(x(1:2,:)-AP3',2) ; ... 
             vecnorm(x(1:2,:)-AP4',2) ]);

%dU = zeros(2,length(Uflow(:,1)));

% % Get the perturbation from solid interactions
dU =  delta1 * ((R + r) > vecnorm(x - app_pos(I,:)')) .* ...
               ((R + r)  - vecnorm(x - app_pos(I,:)')) .* ...
                                  (x - app_pos(I,:)')  ./ vecnorm(x - app_pos(I,:)');

%disp(dU)

% Get the motion resulting from these
dydt = zeros(length(y1),1);
dydt(1:2:end) = Uflow(:,1) + dU(1,:)';
dydt(2:2:end) = Uflow(:,2) + dU(2,:)';

end