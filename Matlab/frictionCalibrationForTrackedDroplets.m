% Hannah Feldstein, Jan Totz, March 2023
% This program plots the loss function between theoretical and measured
% velocity of a droplet as a function of different friction values.

% Load angleVals (output from trackMovieMaker.m), cell matrix containing
% angular position as a function of time.

%%
u0List = 4e-6:4e-6:16e-6; % Velocities to test
muCorrectionList = 0.02:0.02:0.8; % Friction correction factors to test

dropCount = 1;
% Develop color legend
for i = 1:length(u0List)
    colorGrad(i,:) = [(i-1)/length(u0List) 0 (length(u0List)-i)/length(u0List)];
end

dropsToRead = [1 2 3 4  5 6 7 9 10 11 14 15 16 19 20 22 23 24 26 28 31 32];
for dropNum = dropsToRead
    figure
    colorCount = 0;
    thisDrop = angleVals{dropNum};
    for u0 = u0List
        count = 0; % Initialize counter for plotting different speeds
        for muVar = muCorrectionList
            mu = muVar*7.82885e-8;   % Friction coefficient [kg m / s]
            m = 2.20e-13; % in kg
            g = 9.81;   % in m/s^2
            rCM = 2.30e-7; % CoM in m
            
            framerate = 30;         % [1 / s]
            rdroplet = 3.5e-6;         % Radius of the droplet [m]
            fm = 1:length(thisDrop);
            time = fm/framerate; % [s]
            
            % Calculate parameters for each framerate
            theta = thisDrop(:,1);
            phi = thisDrop(:,2);
            
            % Calculate the distance that the bacteria travels, dr
            dr = zeros(1,length(time)-1);
            for i = 1:(length(time)-1)
                dAngle = acos( cos(theta(i))*cos(theta(i+1)) + ...
                sin(theta(i))*sin(theta(i+1))*cos(phi(i+1)-phi(i)) );
                dr(i) = dAngle*rdroplet; % [m]
                dr_velocity(i) = dr(i)/(time(2)-time(1)); % [m / s]
            end
            
            % Calculation of eu, a unit vector that indicates the 
            % direction of the bacterium.
            % Find x,y,z at a given point in time. x and y are given
            % by x_psf and y_psf, respectively:
            x_exp = -rdroplet.*sin(theta).*sin(phi);
            y_exp = rdroplet.*sin(theta).*cos(phi);
            z_exp = rdroplet*cos(theta);

            % For each timestep, find the difference in x,y,z
            dxyz = zeros(length(time)-1,3);
            for i = 1:(length(time)-1)
                dxyz(i,1) = x_exp(i+1)-x_exp(i); % in meters
                dxyz(i,2) = y_exp(i+1)-y_exp(i);
                dxyz(i,3) = z_exp(i+1)-z_exp(i);
            end
            
            % Projection and normalization
            for i = 1:(length(time)-1)
                eu_notnorm(i,:) = sphereProjector(theta(i,:),phi(i,:)) * dxyz(i,:)';
                if norm(eu_notnorm(i,:)) ~= 0
                    eu(i,:) = eu_notnorm(i,:)/norm(eu_notnorm(i,:)); % this line causes trouble, hence the if statement
                else
                    eu(i,:) = eu_notnorm(i,:);
                end
            end
            
            % Calculate Fswim
            muActive = 8.39e-9; % kg m / s
            Fswim = muActive*u0*eu;
            
            % Calculate Fgrav
            a = zeros(length(time)-1,2);
            b = ones(length(time)-1,1);
            F_grav = m*g*rCM/rdroplet*cat(2,a,b);
            
            % Calculate v for a given mu value --> output is a matrix of v for each
            % timestep
            for i = 1:(length(time)-1)
                v(i,:) = (1 / mu) * (sphereProjector(theta(i),phi(i)) * transpose( Fswim(i,:)+F_grav(i,:) ) );
            end
            
            % Calculate the magnitude of velocity
            for i = 1:(length(time)-1)
                v_mag(i) = sqrt( v(i,1)^2 + v(i,2)^2 + v(i,3)^2 );
            end
            
            % Calculate loss function
            lossFunction = 0;
            for i = 1:length(time)-1
                lossFunction = lossFunction + (v_mag(i)-dr_velocity(i)).^2;
            end
    
            % Update counter and save
            count = count + 1;
            loss(count) = lossFunction;
        end
        
        colorCount = colorCount + 1;
        [mini, ind] = min(loss); % Minimum, and index of the minimum value in muVar
        minLoss(colorCount) = muCorrectionList(ind); % The corresponding mu fudge factor for this velocity in question
        semilogy(muCorrectionList,loss,'LineWidth',1,Color=colorGrad(colorCount,:))
        hold on
        xlabel('mu correction factor []','FontSize',16)
        ylabel('loss','FontSize',16)
        title('loss gradient for friction interpolation','FontSize',20)
        
    end
    
    % Form the legend labels
    for i = 1:length(u0List)
        u0ListChar{i} = [num2str(10^6*u0List(i)) , ' microns/s'];
    end
    
    legend(u0ListChar, 'FontSize',14)
    %hold off
    
    for k = 1:length(u0List)
        minLossSave{k}(dropCount) = minLoss(k);
    end
    dropCount = dropCount + 1;
end

function p = sphereProjector(theta,phi)
    p = [cos(theta).^2+sin(theta).^2.*sin(phi).^2, ...
        -cos(phi).*sin(theta).^2.*sin(phi), -cos(theta).*cos(phi).*sin(theta)
    -cos(phi).*sin(theta).^2.*sin(phi), ...
    cos(theta).^2+cos(phi).^2.*sin(theta).^2, -cos(theta).*sin(theta).*sin(phi)
    -cos(theta).*cos(phi).*sin(theta),...
    -cos(theta).*sin(theta).*sin(phi), -cos(theta).*sin(theta).*sin(phi)];
end