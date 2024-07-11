function [Uflow] = surfaceFlow(nTemp1,rot1)

theta = linspace(2*pi/nTemp1,2*pi,nTemp1)';
tangent = [-sin(theta),cos(theta)];
magnitude = zeros(nTemp1,1);

mod = 4/3;
surfF = 300/45;

% if rot1 < pi/2
%     for i = 1:nTemp1
%         if theta(i) >= rot1
%             if theta(i) <= 3/2*pi + rot1
%                 magnitude(i) = -sin(mod*(theta(i)-rot1));
%             else
%                 magnitude(i) = 0;
%             end
%         else
%             magnitude(i) = 0;
%         end
%     end
% else
%     for i = 1:nTemp1
%         if theta(i)>=(rot1-pi/2) && theta(i)<=(rot1)
%             magnitude(i) = 0;
%         else
%             if theta(i)<(rot1-pi/2)
%                 magnitude(i) = -sin(mod*(theta(i)+2*pi-rot1));
%             else
%                 magnitude(i) = -sin(mod*(theta(i)-rot1));
%             end
%         end
%     end
% end

if rot1 <= pi/2
    for i = 1:nTemp1
        if theta(i) >= rot1
            if theta(i) <= 3/2*pi + rot1
                magnitude(i) = -sin(mod*(theta(i)-rot1));
            else
                magnitude(i) = 0;
            end
        else
            magnitude(i) = 0;
        end
    end
else
    for i = 1:nTemp1
        if theta(i)>=(rot1-pi/2) && theta(i)<=(rot1)
            magnitude(i) = 0;
        else
            if theta(i)<=(rot1-pi/2)
                magnitude(i) = -sin(mod*(theta(i)+2*pi-rot1));
            else
                magnitude(i) = -sin(mod*(theta(i)-rot1));
            end
        end
    end
end

Uflow = zeros(nTemp1,2);
Uflow(:,1) = surfF*magnitude.*tangent(:,1);
Uflow(:,2) = surfF*magnitude.*tangent(:,2);

plot(magnitude)
hold on; 
plot(-magnitude(end:-1:1))

%plot( Uflow(:,1).^2 + Uflow(:,2).^2 )

end
