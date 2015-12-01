% This matlab script runs a simulation of a sheet of cells. Individual
% cells divide, die, grow, and secrete morphogens. These actions occur 
% as a result of preprogrammed rules and can use mechanical and chemical
% information as inputs.

clear all
close all

%% Parameter Definitions
% Simulation Parameters
xStep = 0.1;                    % Spatial dimensions of domain [um]
timeStep = 0.01;                % Timestep of simulation [s]
logicStep = 10;                 % Timestep of cell decisions [s]
outputStep = 10;                % Timestep of output [s]
totalTime = 100000;             % Total simulation time [s]
steadyState = 100000;

% Ligand Parameters
spitzCells = [1, 1];            % Number of cells secreting spitz after n divisions [integer]
winglessCells = [1, 1];         % Number of cells secreting wingless after n divisions [integer]
spitzDecayAP = [1e-4, 1e-4];    % Decay coefficient of spitz in A and P compartment [fraction/s]
winglessDecayAP = [1e-4, 1e-4];	% Decay coefficient of wingless in A and P compartment [fraction/s]
spitzDiff = 7e3;                % Diffusivity of spitz [um/s^2]
winglessDiff = 7e3;             % Diffusivity of wingless [um/s^2]
spitzSecretion = 2e-2;          % Spitz secretion over time [AU/s]
winglessSecretion = 2e-2;       % Wingless secretion over time [AU/s]
activeSpiDecay = 1e-2;          % How much does wingless cause spitz to decay [fraction/s]
activeWgDecay = 1e-2;           % How much does spitz cause wingless to decay [fraction/s]

% Cell Parameters
ACells = 8;                     % Number of cells in the A compartment [integer]
PCells = 4;                     % Number of cells in the P compartment [integer]
initialVolume = 10;             % Initial volume of each cell [um]
divisionTimes = [10000];        % Time until nth division [s]
deathSpeed = 0.2;               % Speed at which cells shrink when dying [um/s]

% Probability of apoptosis as a function of chemical and mechanical properties [fraction/s]
deathProb = @(spitz,wingless,engrailed,size) engrailed*(wingless-spitz)*1e-4;

% Growth rate as a function of chemical and mechanical properties [um/s]
growthRate = @(spitz,wingless,engrailed,size) (10 * (spitz^.5 + wingless^.5) - size)/100;

%% Initialization
% Check for invalid inputs
if length(divisionTimes) ~= (length(spitzCells) - 1)
    error('Need to encode spitz secretion size for each division');
end
if length(divisionTimes) ~= (length(winglessCells) - 1)
    error('Need to encode wingless secretion size for each division');
end
if min([spitzCells winglessCells]) < 1
    error('Number of secreting cells must be positive')
end

% Initialize cell array
cellArray = SimulationCellArray(ACells, PCells, xStep);

% Initialize output and update the ligand arrays
outputFrame = 1;
[spitzArray, winglessArray] = ligandArrays(cellArray);
domainSize = length(spitzArray);
outputWindow = domainSize * 3;

%% Solve for steady state PDE solution


%% Simulation Loop
for t = 0:timeStep:totalTime
    % Cell Logic
    if mod(t,logicStep) == 0 && t > steadyState||t==0
        % Update the CELL ligand arrays
        cellArray = updateConcentrations(cellArray, spitzArray, winglessArray);

        % Cell division
        cellArray = divideCellArray(cellArray, t);
          
        % Cell growth and death
        cellArray = growKillCellArray(cellArray, deathProb, growthRate, deathSpeed, logicStep);
        
        % Identify secretion locations
        [spitzSecretionArray, winglessSecretionArray, PCompartment] = findSecretionSites(cellArray, spitzCells, winglessCells);

        % Update the DOMAIN ligand arrays
        [spitzArray, winglessArray] = ligandArrays(cellArray);
        domainSize = length(spitzArray);
    end
    
    % Generate Ligand Arrays
    
    % PDE solving at every timestep
    % Implement Boundary Conditions in a very inefficient way
    spitzArray = periodicBoundaryCondition (spitzArray);
    winglessArray = periodicBoundaryCondition (winglessArray);
    PCompartment = periodicBoundaryCondition (PCompartment);
    spitzSecretionArray = periodicBoundaryCondition (spitzSecretionArray);
    winglessSecretionArray = periodicBoundaryCondition (winglessSecretionArray);
    
    % Solve PDE using lazy Euler (upgrade to better solver later)
    spitzArray = spitzArray + timeStep * ( ...
        -(~PCompartment .* spitzDecayAP(1) + PCompartment .* spitzDecayAP(2)) .* spitzArray + ... % Decay
        -activeSpiDecay .* spitzArray .* winglessArray + ... % Active Decay
        spitzSecretionArray .* spitzSecretion * xStep + ... % Secretion
        spitzDiff .* del2(spitzArray) * xStep^2 ); ... % Diffusion

    winglessArray = winglessArray + timeStep * ( ...
        -(~PCompartment .* winglessDecayAP(1) + PCompartment .* winglessDecayAP(2)) .* winglessArray + ... % Decay
        -activeWgDecay .* spitzArray .* winglessArray + ... % Active Decay
        winglessSecretionArray .* winglessSecretion * xStep + ... % Secretion
        winglessDiff .* del2(winglessArray) * xStep^2 ); ... % Diffusion

    % Undo boundary conditions in a very inefficient way
    spitzArray([1, end]) = [];
    winglessArray([1, end]) = [];
    PCompartment([1, end]) = [];
    spitzSecretionArray([1, end]) = [];
    winglessSecretionArray([1, end]) = [];
    
    % Output
    if mod(t,outputStep) == 0
        volumes = getVolumes(cellArray);
        
        % Show Gradients
        subaxis(2,1,1, 'Spacing', 0, 'Padding', 0, 'Margin', 0)
        x = xStep:xStep:(outputWindow*xStep-xStep);
        domainSize = length(spitzArray);
        repNumber = uint16(outputWindow/domainSize)+1;
        if mod(repNumber,2) == 0
            repNumber = repNumber + 1;
        end
        
        spitzFull = repmat(spitzArray,[1,repNumber]);
        wgFull = repmat(winglessArray,[1,repNumber]);
        leftFull = uint16((length(spitzFull)-outputWindow)/2)+1;
        rightFull = leftFull+outputWindow-2;
        
        plot(x,spitzFull(leftFull:rightFull),x,wgFull(leftFull:rightFull),'LineWidth',2);
        set(gca, 'ColorOrder', [1 0 0; 0 0 1], 'NextPlot', 'replacechildren');
        axis([-inf inf 0 5])
        set(gcf,'color','w');
        axis off
        legend('Spitz', 'Wingless');
        legend boxoff
        
        % Make Canvas
        sheetImage = uint8(ones(400, length(spitzArray), 3))*255;
        
        % Draw Engrailed
        sheetImage(:, logical(PCompartment), [1, 3]) = 0;
        
        % Draw Cell Boundaries
        sheetImage([1:8, 191:200], :, :) = 0;
        sheetImage(:, [1:8], :) = 0;
        positionTemp = 1;
        for i = 1:(length(cellArray.cellArray)-1)
            positionTemp = positionTemp + volumes(i);
            sheetImage(:, uint16((positionTemp-4):(positionTemp+4)), :) = 0;
        end
        
        % Draw Secretion Locations
        sheetImage(201:500,:,:) = 255;
        
        sheetImage(250:350,logical(spitzSecretionArray),2:3) = 0;
        sheetImage(250:350,logical(winglessSecretionArray),1:2) = 0;
        
        sheetImage(250:350,logical(spitzSecretionArray),2:3) = 0;
        sheetImage(250:350,logical(winglessSecretionArray),1:2) = 0;
        
        sheetImage = insertText(sheetImage,[mean(find(spitzSecretionArray)), 300],'Spi','AnchorPoint','Center','TextColor','black','BoxOpacity',0,'FontSize',34);
        sheetImage = insertText(sheetImage,[mean(find(winglessSecretionArray)), 300],'Wg','AnchorPoint','Center','TextColor','black','BoxOpacity',0,'FontSize',34);
        
        
        % Show Image
        subaxis(2,1,2, 'Spacing', 0, 'Padding', 0, 'Margin', 0)
        outputSheet = repmat(sheetImage,[1,repNumber,1]);
        outputSheet = outputSheet(:,(leftFull:rightFull),:);
        outputSheet = insertText(outputSheet,[25, 450],['t = ' datestr(t/3600/24, 'HH:MM:SS')],'AnchorPoint','LeftCenter','TextColor','black','BoxOpacity',0,'FontSize',72);
        outputSheet(:,(end-2):end,:) = 0;
        imshow(outputSheet);
        
        drawnow
        t
        movie(outputFrame) = getframe(gcf);
        outputFrame = outputFrame + 1;
        
    end
end

%% Output Video
movie2avi(movie, 'output.avi', 'compression','none', 'fps', 10);

close all