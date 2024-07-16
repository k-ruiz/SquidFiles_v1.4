% Title: Code to get the fluid flow from a set of Regularized Stokeslets.
% Author: Stephen Williams.

close
tic

%% Add the function files need to run
addpath('functions/')
addpath('classes/')

%% Set parameters
parameters % Set the parameters

N = 50;
U0 = 100/45;
omega = 2*pi/N; 

%% Set channel geometry
stks = getStokesletPositions(rho,system,U0);

%% Solve for the forces
[iS] = getForces(stks,eps_reg);

%% Get the flow over the space

for ii = 1:N

    a = find(stks(:,3) == 8); % Find the Pousielle boundary sections
    Ut = U0*(1+cos(ii*omega))/2; % Set the maximum flow rate in
    stks(a,4:5) = poisuelleFlow(length(a),Ut); % Re-set the Pousielle boundary sections

    a = find(stks(:,3) == 9); % Find the Pousielle boundary sections
    Ut = U0*(1+cos(ii*omega))/2; % Set the maximum flow rate in
    stks(a,4:5) = poisuelleFlow(length(a),Ut); % Re-set the Pousielle boundary sections

    [Uflowx,Uflowy,Uback,omega1] = calculateFlowGrid(stks,iS,x,y,eps_reg);

    hold off;
    Umag = sqrt(Uflowx.^2 + Uflowy.^2);
    imagesc(y,x,Umag); hold on
    c=colorbar;
    c.Limits=[0 10]; % the range that I want
    scatter(stks(:,2),stks(:,1),0.5,'r');
    a = find(stks(:,3) == 2);
    quiver(stks(a(1:4:end),2),stks(a(1:4:end),1),stks(a(1:4:end),5),stks(a(1:4:end),4),'off','k')
    %quiver(stks(x(1:4:end),2),stks(y(1:4:end),1),stks(Uflowx(1:4:end),5),stks(Uflowy(1:4:end),4),'off','k')


    % set(gca,'ylim',[-10 10])
    % set(gca,'xlim',[-12 22])
    axis equal
    saveas(gcf,['outputs/breathCycle/breathCycle_' num2str(ii) '.png'])
    save(['outputs/breathCycle/breathCycle_Uback_omega_' num2str(ii)],'Uback','omega1');
    pause(0.5);
end

%% Make the png into an mp4 file

video = true;

if video == true % AHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    %if havinng problems quality, ffmpeg to convert to better quality. also use handbrake

    v = VideoWriter('breathCycle','MPEG-4');
    v.Quality = 100;
    v.FrameRate = 500;

    open(v);

    for ii = 1:50

        img = imread(['outputs/breathCycle/breathCycle_' num2str(ii) '.png']);
        writeVideo(v,img)

    end

    close(v);

end

toc
