function [Uflow] = surfaceFlow(nTemp1,rot1)

scaler = 300/45;

theta = linspace(2*pi/nTemp1,2*pi,nTemp1);
tangent = [-sin(theta);cos(theta)]';
magnitude = zeros(nTemp1,1);

rot1 = mod(rot1,2*pi);

if rot1 < pi/2 && rot1 ~= 0 

    a = find(theta > rot1);
    b = theta(a) < rot1+3*pi/2;
    b = find(theta(a) < rot1+3*pi/2); b = a(b);
    magnitude(b) = -sin((4/3)*(theta(b)-rot1));

elseif rot1 > pi/2

    a = find(theta < rot1 - pi/2);
    magnitude(a) = -sin((4/3)*(2*pi+theta(a)-rot1));

    b = find(theta > rot1);
    magnitude(b) = -sin((4/3)*(theta(b)-rot1));

elseif rot1 == pi/2

    a = find(theta > pi/2);
    magnitude(a) = -sin((4/3)*(theta(a)-rot1));

else

    a = find(theta < 3*pi/2);
    magnitude(a) = -sin((4/3)*(theta(a)));

end

Uflow = zeros(nTemp1,2);
Uflow(:,1) = scaler*magnitude.*tangent(:,1);
Uflow(:,2) = scaler*magnitude.*tangent(:,2);

end

