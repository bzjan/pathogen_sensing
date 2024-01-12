%% Script for propagation of EM waves beyond the droplet
% Note: Script requires package: statistics and machine learning
% Note: This script must be run in the directory where it is located to execute correctly (use Run instead of Run Section)

% script initialization
close all;                  % close previously open figures
clearvars -except B                     % reset workspace

% choose directory to import from
inputDirectoryNames = ["test_dropletTilt_thetaPhi_analyticDirector_wvl_470_DnLC_0.20" "test_dropletTilt_thetaPhi_analyticDirector_wvl_550_DnLC_0.18" "test_dropletTilt_thetaPhi_analyticDirector_wvl_630_DnLC_0.17"];
channelIndex=2; % Choose which color channel to access
inputDirectoryName = inputDirectoryNames(channelIndex);
inputDirectory = fullfile('.', inputDirectoryName);

% list the polar and azimuthal angles
thetavals = linspace(0,pi/2,40);
phivals = linspace(0,pi/4,20);

% propagation distance zPropagation in nanometers
zPropagation = 1500;               % in nm; seems to agree with experiments
NA = 0.25;                      % numerical aperture of microscope (Olympus 10x air objective)
lambdas = [470,550,630];        % nm
lambda = lambdas(channelIndex); % nm 
gamma = 1.0;                    % gamma correction; 1.0 - no change

outputDirectoryName = "prop_" + inputDirectoryName;
outputDirectory = fullfile('.',outputDirectoryName,sprintf('zDistance_%.1f',zPropagation/1000.0));

% program options
saveH5Q = false;              % save h5 files
showImagesLiveQ = false;      % show non-rescaled images during read in
plotIntensitiesQ = false;      % save single intensity images
useGaussianEnvelopeApproximation = false;
useGammaCorrection = false;
saveMosaicImage = true;
sphericalWavePropagationQ = false;   % spherical wave prop (true) or fresnel prop (false)

% create output directory if it does not exist yet
if not(isfolder(outputDirectory))
    mkdir(outputDirectory);
end

% figure placement options on screen
x0 = 1;
y0 = 1;
width = 421;
height = 400;

% get array size
fieldXpol = get_input_field(round(thetavals(1),3),round(phivals(1),3),inputDirectory,lambda);  % take first file to get array dimensions
[nx,ny] = size(fieldXpol);

% array pre-allocation
integratedIntensity = zeros(length(thetavals),length(phivals));       % Initialize the current intensity matrix
allIntensities = zeros(length(thetavals),length(phivals),nx,ny);  % nx,ny needs to be updated manually!

% discretize the x and y values of intensity images
dx = 70; % dx=dy; [dx] = nm/px (7000 nm / 100 px)
x = (-(nx-1)/2:(nx-1)/2)*dx;
y = x;
[xMat, yMat] = meshgrid(x,y);

% propagator in spherical wave approximation
PointSpreadFunFreeSpace = sphereWaveFun(xMat, yMat, zPropagation, lambda, 0, 0, 0);

% microscope effects
magnification = 1.0;
PSFMicroscopeIncoherent = jincfun(xMat,yMat,2*NA/lambda/magnification, 0, 0).^2;
PSFMicroscopeIncoherent(isnan(PSFMicroscopeIncoherent)) = ((2*NA/lambda/magnification)^2*pi)^2;

% prepare plotting without focus stealing
if showImagesLiveQ
    figIntensityFrameHandle = figure(1);
    axIntensityFrameHandle = axes('Parent',figIntensityFrameHandle);
end

% perform calculations
tic
for q = 1:length(thetavals) 
    for j = 1:length(phivals)
        
        % get electric field data
        theta = round(thetavals(q),3);
        phi = round(phivals(j),3);
        fieldXpol = get_input_field(theta,phi,inputDirectory,lambda);  % get input, single polarization
        
        if useGaussianEnvelopeApproximation
            env = 1200;
            gaussEnv = Gaussfun(xMat, yMat, env);
            fields = fieldXpol.*gaussEnv;
        else
            fields = fieldXpol;
        end
        
        % propagation of EM wave through space
        if zPropagation > 0
            if sphericalWavePropagationQ
                PropagatedField = conv2(PointSpreadFunFreeSpace,fields,'same');         % propagation approximation as a spherical wave
            else
                PropagatedField = fresnel_advance(fields, dx, dx, zPropagation, lambda);
            end
        else
            % no propagation
            PropagatedField = fields;
        end
        
        PropagatedIntensity = abs(PropagatedField).^2;
        ImagedPropagatedIntensity = conv2(PSFMicroscopeIncoherent,PropagatedIntensity,'same');
        
        % live view: image of intensity with axes
        if showImagesLiveQ
            ax = axIntensityFrameHandle; %#ok<*UNRCH> 
            imagesc(ImagedPropagatedIntensity,'Parent',ax);
            colormap(ax,'gray')
            colorbar(ax)
            axis(ax,'square')
            set(ax,'YDir','normal')
            titlename = sprintf('Spherical propagation, $\\theta$=%.3f, $\\phi$=%.3f', theta, phi);
            title(ax,titlename,'interpreter','latex')
        end
        
        % save image dataset to file
        if saveH5Q
            filename = sprintf('data_propagatedIntensity_theta=%.3f_phi=%.3f.png', theta, phi);
            outputfile = fullfile(outputDirectory,[filename,'.h5']);
            h5create(outputfile,"/intensity",size(ImagedPropagatedIntensity))
            h5write(outputfile,"/intensity",ImagedPropagatedIntensity)
        end
        
        % keep data for later plotting with identical colorbar across datasets
        allIntensities(q,j,:,:) = ImagedPropagatedIntensity;

    end
end
toc

%% rescale all intensities so that original min/max of data fall on 0,upperBound
% Note: upperBound is necessary instead of 1, so that max outlier does not distort bulk of data (everything becomes black)
lowerBound = 0.0;
upperBoundPercentage = 97; %  (using 100 for lookup map, was initially 90)
upperBound = max(allIntensities,[],'all')/prctile(allIntensities,upperBoundPercentage,'all');
allRescaledIntensities = rescale(allIntensities,lowerBound,upperBound);

figure(10)
subplot(1,4,1)
    histogram(allIntensities)
    title('raw intensities')
subplot(1,4,2)
    histogram(allRescaledIntensities)
    title('rescaled intensities')
    xlim([0 1])
subplot(1,4,3)
    histogram(allRescaledIntensities.^gamma)
    title('rescaled gamma intensities')
    xlim([0 1])
subplot(1,4,4)
    xVals=linspace(0,1,100);
    plot(xVals,xVals.^gamma)
    hold on
    plot(xVals,xVals)

%% gamma correction
tic

% apply gamma correction
if gamma ~= 1.0
    allRescaledIntensities = allRescaledIntensities.^gamma;
end

for q = 1:length(thetavals)
    for j = 1:length(phivals)
        % scalar data for overview of total intensity
        integratedIntensity(q,j) = sum(allRescaledIntensities(q,j,:,:),'all');
    end 
end
toc

%% Plot of scalar integrated intensity as a function of tilt angles theta and phi I(theta,phi)
if length(thetavals)*length(phivals) > 1 
    
    % rescale total intensities
    normIntegratedIntensity = rescale(integratedIntensity);
    [X,Y] = meshgrid(thetavals,phivals);

    % plot 3d overview
    figure(3);
    img3d = surf(X,Y,normIntegratedIntensity','FaceAlpha',0.5);
    rotate3d on
    xlabel('$\theta \; (^{\circ})$','interpreter','latex')
    ylabel('$\phi \; (^{\circ})$','interpreter','latex')
    zlabel('intensity','interpreter','latex')
    title("propagation distance = " + num2str(zPropagation/1000) + " " + "$\mu$" + "m","interpreter","latex");
    xticks(linspace(0,thetavals(end),7))
    xticklabels(0:15:90)
    yticks(linspace(0,phivals(end),4))
    yticklabels(0:15:45)
    cb = colorbar;
    cb.Label.String = 'Intensity';
    pbaspect([2 1 1])
    outputfile = fullfile(outputDirectory,'intensity_overview_image_3d.png');
    saveas(img3d,outputfile)
    
    % plot 2d overview
    figure(4);
    img2d = imagesc(thetavals,phivals,normIntegratedIntensity');
    xlabel('$\theta \; (^{\circ})$','interpreter','latex')
    ylabel('$\phi \; (^{\circ})$','interpreter','latex')
    title("propagation distance = " + num2str(zPropagation/1000) + " " + "$\mu$" + "m","interpreter","latex");
    xticks(linspace(0,thetavals(end),7))
    xticklabels(0:15:90)
    yticks(linspace(0,phivals(end),4))
    yticklabels(0:15:45)
    cb = colorbar;
    cb.Label.String = 'Intensity';
    set(gca,'YDir','normal')
    pbaspect([2 1 1])
    outputfile = fullfile(outputDirectory,'intensity_overview_image_2d.png');
    saveas(img2d,outputfile)

    % save intensity data overview
    outputfile = fullfile(outputDirectory,'integratedIntensity.mat');
    save(outputfile,'integratedIntensity');
    outputfile = fullfile(outputDirectory,'normIntegratedIntensity.mat');
    save(outputfile,'normIntegratedIntensity')
end

%%%%%%% post analysis

%% compute intensity plot range
intensityLimits = [0.0 1.0]; 
% tool for finding the right intensity limits manually
figHisto = figure(5);
histogram(allRescaledIntensities)
xlim([-1 2])
xline(intensityLimits,'-r',{'I_{min}','I_{max}'})

%% create overview mosaic image of all single intensity images
if length(thetavals)*length(phivals) > 1 
    nx = 101;
    ny = 101;
    intensityMosaic = zeros(length(thetavals)*nx,length(phivals)*ny);
    for q = 1:length(thetavals)
        for j = 1:length(phivals)
            intensityMosaic((1:nx)+nx*(q-1),(1:ny)+ny*(j-1)) = squeeze(allRescaledIntensities(q,j,:,:));
        end
    end
    figIntensitiesMosaic = figure('Position', [10 10 1800*2 900*2]);
    ax = axes('Parent',figIntensitiesMosaic);
    img = imagesc(intensityMosaic');
    colormap('gray')
    pbaspect([2 1 1])
    set(ax,'YDir','normal')
    xlabel('$\theta \; (^{\circ})$','interpreter','latex')
    ylabel('$\phi \; (^{\circ})$','interpreter','latex')
    title("propagation distance = " + num2str(zPropagation/1000) + " " + "$\mu$" + "m","interpreter","latex");
    xticks(linspace(1,length(thetavals)*nx,7))
    xticklabels(0:15:90)
    yticks(linspace(1,length(phivals)*ny,4))
    yticklabels(0:15:45)
    clim(ax,intensityLimits)
    outputfile = fullfile(outputDirectory,'all_intensity_images_overview.png');
    if saveMosaicImage
        exportgraphics(ax,outputfile,'Resolution',1200)
    end
end

if plotIntensitiesQ
    tic
    
    % plot images
    % get figure and axes handles for plotting without focus stealing
    figIntensityFrameHandle = figure(1);
    clf(figIntensityFrameHandle,'reset');
    axIntensityFrameHandle = axes('Parent',figIntensityFrameHandle);
    
    figIntensityNoFrameHandle = figure(2);
    clf(figIntensityNoFrameHandle,'reset');
    axIntensityNoFrameHandle = axes('Parent',figIntensityNoFrameHandle);
    
    % iterate over all datasets
    for q = 1:length(thetavals) 
        for j = 1:length(phivals)
            theta = round(thetavals(q),3);
            phi = round(phivals(j),3);
            ImagedPropagatedIntensity = squeeze(allRescaledIntensities(q,j,:,:));
            
            % image without axes
            ax = axIntensityNoFrameHandle;
            imgNoAxes = imagesc(ImagedPropagatedIntensity,'Parent',ax);
            colormap(ax,'gray')
            clim(ax,intensityLimits)
            axis(ax,'square')
            set(ax,'YDir','normal')
            set(ax, 'visible', 'off')
            filename = sprintf('plotNoAxes_propagatedIntensity_%.0f_theta=%.3f_phi=%.3f.png',lambda, theta, phi);
            outputfile = fullfile(outputDirectory,filename);
            exportgraphics(ax,outputfile)
        end
    end
    toc
end

function fieldXpol = get_input_field(theta,phi,inputDirectory,wvl)
    filename = sprintf('afterDropletEField_%d_theta_%.3f_phi_%.3f.h5',wvl,theta,phi); % Original file naming system
    pthfn = fullfile(inputDirectory,filename);              % pth and filename of data
    eFields = h5read(pthfn,'/Dataset1');                    % contains Ex,Ey,x,y (2,2,100,100)
    eFields = permute(eFields,[4,3,2,1]);                   % h5 import flips array order of mathematica
    fieldXpol = eFields(:,:,2,1) + 1i*eFields(:,:,2,2);     % output field of single polarization to propagate
end

% Spherical wave kernel
function output = sphereWaveFun(x, y, z, lambda, xoffset, yoffset, zoffset)
    arg = 1i*2*pi/lambda * sqrt((x-xoffset).^2 + (y-yoffset).^2 + (z-zoffset).^2);
    output = exp(arg)./arg;
end

% Jinc for microscope blurring kernel
function output = jincfun(x, y, r0, offsetx, offsety)
    rxy = sqrt((x-offsetx).^2 + (y-offsety).^2);
    output = r0*besselj(1,2*pi*r0*rxy)./rxy;
end

% Gauss to approximate lens effect of droplet
function output = Gaussfun(x,y, width)  %#ok<DEFNU>
    output = exp(-(x.^2 + y.^2)/width^2);
end

function eProp = fresnel_advance(U0, dx, dy, z, lambda)
	% The function receives a field U0 at wavelength lambda
	% and returns the field U after distance z, using the Fresnel
	% approximation. dx, dy, are spatial resolution.

	k=2*pi/lambda;
	[ny, nx] = size(U0); 

	Lx = dx * nx;
	Ly = dy * ny;

	dfx = 1./Lx;
	dfy = 1./Ly;

	u = ones(nx,1)*((1:nx)-nx/2)*dfx;    
	v = ((1:ny)-ny/2)'*ones(1,ny)*dfy;   

	O = fftshift(fft2(U0));

	H = exp(1i*k*z).*exp(-1i*pi*lambda*z*(u.^2+v.^2));  
	eProp = ifft2(ifftshift(O.*H));
end

