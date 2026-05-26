%% View Timetraces

clear,close all
load('angleVals.mat') % run trackDroplets.m + save variable to get angleVals.mat

frameRate = 30; 
viewAngles = 1;

% find datasets with highest thetas in decreasing order; store maximum
% values in maxVal and array indices in maxIndex
maxAngleVals = zeros(length(angleVals),1);
for ii = 1:length(angleVals)
    maxAngleVals(ii) = max(angleVals{ii}(:,1));
end
for ii = 1:10 % set this value higher if needed, limited at 38
    [maxVal(ii), maxIndex(ii)] = max(maxAngleVals);
    maxAngleVals(maxIndex(ii)) = 0; 
end

ii=8; % variable that selects the dataset number from maxIndex (1-10)

t = (0:1:size(angleVals{maxIndex(ii)},1)-1)/frameRate; % time array 

% apply median filter to noisy data 
filterOrder = 15; % median filter order 
filteredTheta = medfilt1(angleVals{maxIndex(ii)}(:,1),filterOrder);
unwrappedPhi = unwrap(angleVals{maxIndex(ii)}(:,2));
filteredPhi = medfilt1(unwrap(angleVals{maxIndex(ii)}(:,2)),filterOrder);

% smooth median-filtered data 
sfactor=15;
filteredTheta2 = smooth(filteredTheta,sfactor);
filteredPhi2 = smooth(filteredPhi,sfactor);

dt = t(2)-t(1); % dt

% time derivatives of theta and phi for unfiltered, medianfiltered, and
% medianfiltered + smoothed data
dThetadtUnfilt = 1/dt * (angleVals{maxIndex(ii)}(2:end,1) - angleVals{maxIndex(ii)}(1:end-1,1)); 
dThetadtFilt = 1/dt * (filteredTheta(2:end) - filteredTheta(1:end-1));
dThetadtFilt2 = 1/dt * (filteredTheta2(2:end) - filteredTheta2(1:end-1));

dPhidtUnfilt = 1/dt * (angleVals{maxIndex(ii)}(2:end,2) - angleVals{maxIndex(ii)}(1:end-1,2)); 
dPhidtFilt = 1/dt * (filteredPhi(2:end) - filteredPhi(1:end-1));
dPhidtFilt2 = 1/dt * (filteredPhi2(2:end) - filteredPhi2(1:end-1));

lLim = 0; % lower time axis limit for plotting 
uLim = max(t);  % high time axis limit for plotting 

if viewAngles
    figure('Units','normalized', 'Position',[0 0 0.5 1])
        tiledlayout(4,1)
        nexttile(1)
            plot(t, angleVals{maxIndex(ii)}(:,1));
            hold on
            plot(t, filteredTheta,'LineWidth', 1);
            plot(t, filteredTheta2,'LineWidth', 1.5);
            legend('raw','median','smooth')
            xlim([lLim uLim])
            ylabel('\theta [rad]')
        nexttile(2)
            plot(t(1:end-1), dThetadtUnfilt);
            hold on 
            plot(t(1:end-1), dThetadtFilt, 'LineWidth', 1);
            plot(t(1:end-1), dThetadtFilt2, 'LineWidth', 1.5);
            xlim([lLim uLim])
            ylabel('d\theta/dt [rad/s]')
        nexttile(3)
            plot(t, unwrappedPhi);
            hold on
            plot(t, filteredPhi, 'LineWidth', 1);
            plot(t, filteredPhi2,'LineWidth', 1.5);
            xlim([lLim uLim])
            ylabel('\phi [rad]')
        nexttile(4)
            plot(t(1:end-1), dPhidtUnfilt);
            hold on 
            plot(t(1:end-1), dPhidtFilt, 'LineWidth', 1);
            plot(t(1:end-1), dPhidtFilt2, 'LineWidth', 1.5);
            xlim([lLim uLim])
            ylabel('d\phi/dt [rad/s]')
            xlabel('t [s]')
end

%% Loop over all specified values - obtain mean and standard deviation
close all

% Load parameters
rdroplet = 3.5e-6;  
m = 2.20e-13; % in kg
g = 9.81;   % in m/s^2
rCM = 2.30e-7; % CoM in m
muActive = 8.39e-9; % kg / s

% loop over nx3 matrix, where columns are (ii,tStart,tEnd)
params = [2 6.5 7.5
    2 60.5 61.23
    3 22 23.3
    3 33.9 35
    4 1.0 2
    5 18.1 19.0
    5 31.1 31.8
    8 3.16 4.1
    8 8.46 9.15
    1 6.3 7.6
    5 9.5 10.3];

muFactors = 0.025:0.025:1;

for i = 1:length(params)
    % Settings
    clear mu

    ii = params(i,1); % select droplet
    t1 = params(i,2); % start time (s)
    t2 = params(i,3); % end time (s)
    
    % Extract data
    thetaVals = angleVals{maxIndex(ii)}(:,1); % time array 
    phiVals = angleVals{maxIndex(ii)}(:,2); % time array 
    fms = round(t1*frameRate):round(t2*frameRate);
    thetaTimeseries = thetaVals(fms); % get theta values over chosen frames
    phiTimeseries = phiVals(fms);
    
    % Smooth data
    % apply median filter to noisy data 
    filterOrder = 10; % median filter order - 10
    filteredTheta = medfilt1(thetaTimeseries,filterOrder);
    filteredPhi = medfilt1(phiTimeseries,filterOrder);
    % smooth median-filtered data 
    sfactor=15; % # of values to avg over - 15
    filteredTheta2 = smooth(filteredTheta,sfactor);
    filteredPhi2 = smooth(filteredPhi,sfactor);

    % Calculate theta at max value/leveling off
    thetaMax = max(thetaTimeseries);
    F_grav = m*g*rCM/(rdroplet.^2)*sin(thetaMax);
    u = F_grav*rdroplet/muActive; %u0 if torque due to bacterium swimming is balance by droplet torque
    uVals(i) = u;

    figure
    subplot(3,1,1)
        hold on
        box on
        plot(fms/frameRate,filteredTheta2*180/pi,'LineWidth',2)
        plot(fms/frameRate,thetaTimeseries*180/pi,'LineWidth',2)
        xlabel('time (s)')
        ylabel('theta (degrees)')
        if i == 8
            xlim([3.2 4.1])
        end

    for j = 1:length(muFactors)
        % Solve differential equationsyms y(t) a
        correctionFactor = muFactors(j);
        muInput = correctionFactor*7.82885e-8;

        timeAxis = t1:(1/30):t2;
        tspan = timeAxis-t1;
        y0 = filteredTheta2(1);
        [t,y] = ode45(@(t,y) (muActive*u/rdroplet - m*g*rCM/rdroplet^2*sin(y))/muInput, tspan, y0);

        if length(y)>length(filteredTheta2)
            y = y(1:(end-1));
        elseif length(y)<length(filteredTheta2)
            filteredTheta2 = filteredTheta2(1:(end-1));
        end

        totalLoss = (y-filteredTheta2).^2;
        lossCutoff = 60; % how long to read signal for (# of frames)
        if length(totalLoss) > lossCutoff
            loss(j) = sum(totalLoss(1:lossCutoff)); % take first half second of data
        else
            loss(j) = sum(totalLoss(1:end-3));
        end

    end

    subplot(3,1,2)
        hold on
        box on
        plot(fms/frameRate,phiTimeseries*180/pi,'LineWidth',1)
        xlabel('time (s)')
        ylabel('phi (degrees)')
        ylim([0 360])
    subplot(3,1,3)
        box on
        plot(muFactors,loss/max(loss))
        xlabel('correction factor')
        ylabel('normalized loss [degrees^2]')

    % Find the best correction factor
    [M,indexMu] = min(loss);
    bestFitCorrectionFactors(i) = muFactors(indexMu);
    
    
    muInput = bestFitCorrectionFactors(i)*7.82885e-8;
    [t,y] = ode45(@(t,y) (muActive*u/rdroplet - m*g*rCM/rdroplet^2*sin(y))/muInput, tspan, y0);
    subplot(3,1,1)
        plot(timeAxis,y*180/pi,'LineWidth',1)
        box on
        legend('filtered','unfiltered','best fit')
    
    
end

disp('mean correction factor:')
disp(mean(bestFitCorrectionFactors))
disp('standard dev:')
disp(std(bestFitCorrectionFactors))


