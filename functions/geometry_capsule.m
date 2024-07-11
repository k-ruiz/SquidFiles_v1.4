% function to initialize the outer boundary / mantle
% walls A, B, C should have stks(:,3) = 8
% walls D, E should have stks(:,3) = 9

% walls A and C are symmetrical to each other, side walls
% walls D and E are symmetrical to each other, bottom walls
% wall B is the top wall

function [caps_stks] = geometry_capsule(rho2,Pty2,Ltot2,side2,bott2,top2)

    %% construction of the top boundary
    phi = linspace(0,pi,ceil(top2*rho2*pi)); % gets the angles to parameterise the surface of the top boundary

    wallB = [cos(phi);sin(phi);zeros(ceil(top2*rho2*pi),1)']'; % gets the coordinates of the half-circle
    wallB = floor(top2)*wallB + [0,side2,0]; % shifts the half-circle to the top of the squid

    %% construction of the bottom boundaries
    wallD = [linspace(bott2,top2-1,floor(top2*rho2));ceil(Pty2-Ltot2)*ones(floor(top2*rho2),1)';8*ones(floor(top2*rho2),1)']';
    wallD = wallD(2:end-1,:); % removes the overlapping points with the other walls

    wallE = [linspace(-bott2,-top2+1,floor(top2*rho2));ceil(Pty2-Ltot2)*ones(floor(top2*rho2),1)';9*ones(floor(top2*rho2),1)']';
    wallE = wallE(2:end-1,:); % removes the ovverlapping points with the other walls

    %% construction of the side boundaries 
    wallA = [floor(top2)*ones(ceil(2*side2*rho2),1)';linspace(Pty2-Ltot2,side2,ceil(2*side2*rho2));zeros(ceil(2*side2*rho2),1)']';
    wallA = wallA(1:end-1,:); % removes the overlapping points with the other walls

    wallC = wallA;
    wallC(:,1) = -wallC(:,1); % reflect acorss the x axis

    %% combining all the arrays
    caps_stks = [wallA;wallB;wallC;wallD;wallE]; % combines all the boundaries

end
