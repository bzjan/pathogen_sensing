% Hannah Feldstein, Jan Totz
% This script plots the loss function between theoretical and measured
% polar angle velocity of a droplet as a function of different friction 
% values.

% Load angle values
load('angleVals.mat')

% Choose certain plots (1=yes,0=no)
plotIndividual=1; % Plot individual loss functions
plotDropVelocity=0; % Plot droplet velocity
smoothAngles = 0; % Apply a smoothing function to angles
plotAngles=0; % Plot angles for each droplet
sfactor = 5; % Smoothing factor, default value in MATLAB is 5
figC=0;

u0List = 6e-6:4e-6:18e-6; % Swimming speed values to test [m/s]
muCorrectionList = 0.02:0.02:4; % Friction correction factors to test

% Develop color legend
colorCount = 1; % Counter for plot color
dropCount = 1;
for i = 1:length(u0List)
    colorGrad(i,:) = [(i-1)/length(u0List) 0 (length(u0List)-i)/length(u0List)];
end

% Only read droplets with a long enough time signal
dropsToRead = [];
ii = 1;
for i = 1:length(angleVals)
    if length(angleVals{i}) > 2000
        dropsToRead(ii) = i;
        ii=ii+1;
    end
end

for dropNum = dropsToRead
    fprintf('dropNum = %.2e\n\n',dropNum)
    % figure
    colorCount = 0;
    % Run this code on each of the droplets in question, i.e. through all cells
    % of the angleVals cell matrix
    thisDrop = angleVals{dropNum};
    for u0 = u0List
        count = 0; % initialize counter for plotting different speeds
        clear F_grav FswimTheta FswimPhi vTheta vPhi vThetaTheory vPhiTheory
        for muTest = muCorrectionList
            % Testing parameters
            mu = muTest*7.82885e-8;   % friction coefficient in kg / s
            m = 2.20e-13; % in kg
            g = 9.81;   % in m/s^2
            rCM = 2.30e-7; % CoM in m

            % Main script
            framerate = 30;         % in Hz
            rdroplet = 3.5e-6;         % Radius of the droplet in m
            fm = 1:length(thisDrop);
            time = fm/framerate; % Time in seconds
            tTot = length(time)-1;
            
            % Calculate parameters for each framerate
            thetaRaw = thisDrop(:,1);
            phiRaw = thisDrop(:,2);
            phiRaw = unwrap(phiRaw);
            if smoothAngles
                theta=smooth(thetaRaw,sfactor);
                phi=smooth(phiRaw,sfactor);
            else
                theta=thetaRaw;
                phi=phiRaw;
            end

            % Calculate Fgrav
            F_grav = transpose(m*g*rCM/(rdroplet.^2)*sin(theta));

            % Calculate Fswim
            muActive = 8.39e-9; % kg m / s
            vCounterFg = 0;
            for i=1:tTot
                % Calculate the change in polar and azimuthal angles
                dTheta(i) = theta(i+1)-theta(i); dPhi(i) = phi(i+1)-phi(i);
    
                % Calculate the angular velocities
                vTheta(i) = dTheta(i)/(time(2)-time(1)); vPhi(i) = dPhi(i)/(time(2)-time(1));

                if vPhi(i) == 0 % If both are 0, the bacteria is still swimming in the theta direction to counter Fg
                    FswimTheta(i) = muActive*u0/rdroplet;
                    FswimPhi(i) = 0;
                else % case where neither = 0
                    vCounterFg = F_grav(i)/muActive; % swimming velocity that counters F_grav
                    omega = atan(vPhi(i)/(vTheta(i)+vCounterFg)); % angle that describes relative theta and phi velocities
                    FswimTheta(i) = muActive*u0/rdroplet*cos(omega);
                    FswimPhi(i) = muActive*u0/rdroplet*sin(omega);
                end

                % Correct the sign to correspond to direction
                FswimTheta(i)=abs(FswimTheta(i))*sign(vTheta(i)+vCounterFg);
                FswimPhi(i)=abs(FswimPhi(i))*sign(vPhi(i));
            end

            
            % Calculate v for a given friction (mu) value
            % Calculate v for this given timestep
            vThetaTheory = (1/mu) * (FswimTheta - F_grav(1:(end-1)));
            vPhiTheory = (1/mu) * FswimPhi;
            
            % Calculate loss function
            lossFunctionTheta = 0;
            lossFunctionPhi = 0;
            for i = 1:tTot
                lossFunctionTheta = lossFunctionTheta + (vTheta(i)-vThetaTheory(i)).^2;
                lossFunctionPhi = lossFunctionPhi + (vPhi(i)-vPhiTheory(i)).^2;
            end
    
            % Update counter and save
            count = count + 1;
            lossTheta(count) = lossFunctionTheta;
            lossPhi(count) = lossFunctionPhi;
        end

            colorCount = colorCount + 1;
            [miniTheta, indTheta] = min(lossTheta); % Minimum, and index of the minimum value in mu list
            [miniPhi, indPhi] = min(lossPhi); % Minimum, and index of the minimum value in mu list
            minLossTheta(colorCount) = muCorrectionList(indTheta); % The corresponding mu correction factor for this velocity in question
            minLossPhi(colorCount) = muCorrectionList(indPhi);
        if plotIndividual
            if u0==min(u0List)
                figure
            end
            % semilogy(muCorrectionList,lossTheta,'LineWidth',1,Color=colorGrad(colorCount,:))
            semilogy(muCorrectionList,lossTheta,'LineWidth',2)
            hold on
            xlabel('mu correction factor []','FontSize',16)
            ylabel('loss','FontSize',16)
            title('loss gradient for theta friction interpolation','FontSize',20)
        end
    end
    
    % Form the legend labels
    for i = 1:length(u0List)
        u0ListChar{i} = [num2str(10^6*u0List(i)) , ' microns/s'];
    end
    legend(u0ListChar, 'FontSize',14)

    for k = 1:length(u0List)
        minLossSaveTheta{k}(dropCount) = minLossTheta(k);
        minLossSavePhi{k}(dropCount) = minLossPhi(k);
    end
    dropCount = dropCount + 1;


% Plot the x,y,z, components of drop velocity as a function of time
if plotDropVelocity
    figure
    subplot(3,1,1)
        plot(time(1:end-1),v(:,1))
    subplot(3,1,2)
        plot(time(1:end-1),v(:,2))
    subplot(3,1,3)
        plot(time(1:end-1),v(:,3))
end

% Plot the theta and phi to observe the effect of smoothing

if plotAngles
    figC=figC+1;
    figure(figC)
    
    subplot(2,1,1)
        hold on
        plot(time,thetaRaw,'r')
        plot(time,theta,'b')
        xlabel('time (s)')
        ylabel('polar angle (rads)')
        hold off
    subplot(2,1,2)
        hold on
        plot(time,phiRaw,'r')
        plot(time,phi,'b')
        xlabel('time (s)')
        ylabel('azimuthal angle (rads)')
        hold off
end


end



%% plot individual curves theta
A = minLossSaveTheta{1};
B = minLossSaveTheta{2};
C = minLossSaveTheta{3};
D = minLossSaveTheta{4};

n = 10;

h = histfit(A,n,'kernel');
set(gcf,'Visible','off')
CurveXA = h(2).XData;
CurveYA = h(2).YData;

h = histfit(B,n,'kernel');
set(gcf,'Visible','off')
CurveXB = h(2).XData;
CurveYB = h(2).YData;

h = histfit(C,n,'kernel');
set(gcf,'Visible','off')
CurveXC = h(2).XData;
CurveYC = h(2).YData;

h = histfit(D,n,'kernel');
set(gcf,'Visible','off')
CurveXD = h(2).XData;
CurveYD = h(2).YData;

figure(61)
hold on
plot(CurveXA,CurveYA,'LineWidth',2)
plot(CurveXB,CurveYB,'LineWidth',2)
plot(CurveXC,CurveYC,'LineWidth',2)
plot(CurveXD,CurveYD,'LineWidth',2)
legend(u0ListChar,'FontSize',16)
xlabel('correction factor','FontSize',12)
ylabel('occurance','FontSize',12)
title('Polar (theta)','FontSize',18)
box on
hold off
