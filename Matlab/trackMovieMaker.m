% Hannah Feldstein, Jan Totz, March 2023
% Requirements: Image Processing Toolbox
% This program is meant to automatically track droplets, and to
% calculate the distance between the center of the droplet and the
% bright spot that appears in linearly polarized light.

% 0) Read all frames from movie
% 1) Find the location (center) of all droplets and build trajectories
% 2) Determine where on the droplet the bright spot is located and
% calculate the orientation (polar and azimuthal angle) of the droplet

% options
debugQ = true; % Runs through selected range of frames
% debugQ = false; % Runs through all frames
makeMovie = true; % Saves a video
makeMovie = false;

if makeMovie
    % Select which droplets to record in the video - each droplet is
    % automatically assigned a number
    dropletsToFilm = [3 17];
end


%% 0) Read video
tic
videoStream = VideoReader("Supplementary Video 4.mp4"); % Read the video
nFrames = videoStream.NumFrames;
% frameRate = videoStream.FrameRate;
frameRate = 30; % Specify framerate if video was not taken at 1x

if debugQ
    % Select range of frames to analyze, if debugQ
    startFrame = 300;
    endFrame = 500;
    frameRange = [startFrame,endFrame];
    allFrames = read(videoStream,frameRange);
else
    startFrame = 1;
    endFrame = nFrames;
    allFrames = read(videoStream);
end
[ny,nx,~,numFrames] = size(allFrames);

time = toc;
fprintf("Time to read frames: %.1f s\n",time);

%% 1) Find droplets
tic
showCircles = false; % shows locations of circles for debugging
radiusMinMax = [7,20];    % in px; find via drawline()
for t = 1:numFrames
    [centers,radii] = imfindcircles(allFrames(:,:,:,t),radiusMinMax, ...
        "ObjectPolarity","dark", ...
        "Sensitivity",0.89, ...
        "EdgeThreshold",0.1, ...
        "Method","twostage");
    allCenters{t} = centers; % Center pixel of each found circle

    if showCircles
        figure(1);
        imshow(allFrames(:,:,:,t))
        viscircles(centers,radii)
    end
end

time = toc;
fprintf("Time to find droplets: %.1f s\n",time);

%% Build trajectories from droplet positions

dropTraj = {}; % Cellular array
% If the nearest next-frame droplet is more than this distance from the original droplet
% in the next frame, the trajectory will be cut short
thresholdDist = 10;

[numDrops,xyCoords] = size(allCenters{1});
for i = 1:numDrops
    dropTraj{i}(1,:) = allCenters{1}(i,:); % Each cell corresponds to a droplet

    % Cycle through all of the found droplets
    for j = 1:(endFrame-startFrame)
        closestPoint = nearestDrop(dropTraj{i}(j,:),allCenters{j+1}); % closest point
        if dist(dropTraj{i}(j,1),dropTraj{i}(j,2),closestPoint(1),closestPoint(2)) <= thresholdDist
            dropTraj{i}(j+1,:) = closestPoint;
        else
            break
        end

    end
end

%% Plot the trajectories of the droplets
plotTraj = true;
figure(2)
hold on
if plotTraj
    for i = 1:length(dropTraj)
        scatter(dropTraj{i}(:,1),dropTraj{i}(:,2))
    end
end
hold off

%% 2) Find bright spot and orientation of each droplet

% Set radius of the mask - adjust based on droplet size in video
rEye = 10;
tic
% Save a cellular matrix containing values of theta and phi (to be
% referenced later)
thetaPhiMasks = {}; % Matrix containing theta and phi values
thetaRes = 90; % Number of theta values
phiRes = 40; % Number of phi values
thetaVals = linspace(0,pi/2,thetaRes);
phiVals = linspace(0,2*pi,phiRes);
for i = 1:thetaRes
    theta = thetaVals(i);
    for j = 1:phiRes
        phi = phiVals(j);
        % Create a mask
        thetaPhiMasks{i,j} = generateMask(theta,phi,rEye);
    end
end
time = toc;
fprintf("Time to generate masks: %.1f s\n",time);

% Initialize a cell matrix that contains [theta,phi] for each droplet at each timestep
angleVals = {};

% Add a zero padding factors
zp = 5*rEye;

% Set up the movie
if makeMovie
    v = VideoWriter('newfile.avi');
    v.FrameRate = 30;
    open(v)
    h = figure;
end

tic

% Follow each individual droplet
for dropNum = 1:length(dropTraj)
    disp(dropNum)
    % Look at each timestep
    %disp(dropNum/length(dropTraj))
    drop = dropTraj{dropNum}; % Nx2 matrix containing x-y coords
    [h,w] = size(drop);
    for td = 1:h
        % For each point in time, find the droplet and take the convolution
        % with different theta and phi
        % Read the frame
        imFrameRaw = allFrames(:,:,:,td);

        % Zero pad the frame
        imFrame = padarray(imFrameRaw,[zp zp],1,'both');

        yTest = round(drop(td,1))+zp; % Add padding
        xTest = round(drop(td,2))+zp;
        imDrop = imFrame((xTest-rEye):(xTest+rEye),(yTest-rEye):(yTest+rEye));
        imDropBlurred = imgaussfilt(imDrop,1.5); % Blur to denoise
        imDropBlurredBWNotPadded = im2bw(imDropBlurred,0.18); % Binarize

        % Cycle through all values of theta/phi, and take
        % the convolution. Whichever theta/phi correspond to the highest
        % avg are saved as theta/phi values of the droplet.
        bestTheta = 0; % Initialize the theta/phi corresponding to best match
        bestPhi = 0;
        bestConvVal = 0;
        for i = 1:thetaRes
            for j = 1:phiRes
                % Take the convolution
                mask = thetaPhiMasks{i,j};
                % Zero pad the mask by a factor of 7
                imDropBlurredBW = padarray(imDropBlurredBWNotPadded,[7 7],1,'both');
                convVal = mean(mean(mask.*imDropBlurredBW));
                
                if convVal > bestConvVal
                    bestConvVal = convVal;
                    bestTheta = thetaVals(i);
                    bestPhi = phiVals(j);
                end
            end
        end
        % Add best theta and phi values to the matrix
        angleVals{dropNum}(td,1) = bestTheta;
        angleVals{dropNum}(td,2) = bestPhi;
        
        if makeMovie && ismember(dropNum,dropletsToFilm)
            titlestring = {strcat('from video, i = ',num2str(dropNum)), strcat('t = ',num2str(round(td/30)), 's') };
            subplot(2,2,1)
                imshow(imDrop)
                title('Imaged Droplet','FontSize',14)
            subplot(2,2,2)
                imshow(generateMask(bestTheta,bestPhi,rEye))
                title('Modeled Droplet','FontSize',14)
            subplot(2,2,3)
                plot((1:length(angleVals{dropNum}(:,1)))/frameRate,180/pi*angleVals{dropNum}(:,1),'LineWidth',0.75)
                axis([0, (endFrame-startFrame)/frameRate, 0, 90])
                title('Polar Angle','FontSize',14)
                axis square
                xlabel('Time (s)','FontSize',14)
                ylabel('$\theta^{o}$','fontsize',14,'interpreter','latex')
            subplot(2,2,4)
                plot((1:length(angleVals{dropNum}(:,2)))/frameRate,180/pi*angleVals{dropNum}(:,2),'LineWidth',0.75)
                axis([0, (endFrame-startFrame)/frameRate, 0, 360])
                axis square
                title('Azimuthal Angle','FontSize',14)
                xlabel('Time (s)','FontSize',14)
                ylabel('$\phi^{o}$ ','FontSize',14,'interpreter','latex')
            F = getframe(h);
            writeVideo(v,F.cdata)
        end

    end
end

if makeMovie
    close(v);
end

time = toc;
fprintf("Time to track droplet orientation: %.1f s\n",time);

function d = dist(x0,y0,x1,y1)
% A simple distance function, neglecting square root
    d = (x1-x0)^2 + (y1-y0)^2;
end

function closest = nearestDrop(originalDrop,nextFrameDroplets)
% originalDrop: [x0, y0], coordinates of the original droplet
% nextFrameDroplets: Nx2, coordinates of all the other identified droplets
% in next frame
    [numRow,numCol] = size(nextFrameDroplets);
    minDist = dist(originalDrop(1),originalDrop(2),nextFrameDroplets(1,1),nextFrameDroplets(1,2)); % minimum distance
    besti = 1; % initialize the index corresponding to the closest droplet in the next frame
    for i = 2:numRow
        d = dist(originalDrop(1),originalDrop(2),nextFrameDroplets(i,1),nextFrameDroplets(i,2)); % find distance
        if d < minDist
            minDist = d; % Reset minimum distance
            besti = i;
        end
    end
    closest = nextFrameDroplets(besti,:);
end

function mask = generateMask(theta,phi,rDrop,rBright)
% This function generates a binary mask, based on the angular position
% of the droplet.
% Inputs:
% thetha - polar angle, radians
% phi - azimuthal angle, radians
% rDrop - radius of the mask of the droplet
% rBright - radius of the mask of the bright spot
% Outputs:
% mask - mask of the droplet based on corresponding theta and phi [NxN]
    
    if nargin < 4
        rBright = 5;
    end
    if nargin < 3
        rDrop = 10;
    end

    % Create the mask
    sizeMask = round(rDrop*3.5);
    maskEdge = -floor(sizeMask/2):1:floor(sizeMask/2);
    mask = (sqrt(maskEdge.^2 + maskEdge'.^2) >= rDrop);

    % Form bright spot corresponding to the angular position
    distFromCenter = rDrop*sin(theta);
    xp = round(sizeMask/2)+round(distFromCenter*cos(phi)); % x-coord of the bright center
    yp = round(sizeMask/2)-round(distFromCenter*sin(phi)); % y-coord of the bright center
    for i = -rBright:rBright
        xPos = xp+i;
        for j = -rBright:rBright
            yPos = yp+j;
            if (xPos > 0) && (yPos > 0)
                if sqrt((xPos-xp)^2+(yPos-yp)^2) <= rBright
                    mask(xPos,yPos) = 1;
                end
            end
        end
    end
end