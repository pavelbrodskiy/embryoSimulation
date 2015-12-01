%% Definitions
% Simulation Parameters
domainSize = 200;       %microns
xStep = 0.1;            %microns
totalTime = 1000;       %s
timeStep = 0.001;       %s
logicStep = 0.1;        %s
outputTime = 0.1;         %s

% Spitz Parameters
spitzSecretion = 1;     %uM
spitzDecayA = 0.05;      %Percent decayed per second
spitzDecayP = 0.05;      %Percent decayed per second
spitzDiff = 222;        %uM/s^2

% Cell Parameters
cellSize = 5;
optimalCellSize = 2;

deathProb = @(spitz,size) (0.5-spitz)*0.01/cellSize;
divProb = @(spitz,size) max([0,(cellSize-optimalCellSize)/1000]);


%% Initialization

x = 0:xStep:domainSize;
arrayLength=length(x);
blankArray = ones(1,arrayLength);
secretionArray = blankArray*0;
secretionArray(uint16(secretionCoord/xStep)) = 1;
spitzArray = blankArray*0;

cellBoundaries = blankArray*0;
cellBoundaries(1:cellSize/xStep:arrayLength) = 1;
cellBoundaries(end) = 1;

cellNumber = length(find(cellBoundaries));
j = 8;
for i = 1:cellNumber
    j = j - 1;
    engrailed(i) = j > 4;
    if j == 0
        j = 8;
    end
end

%% Simulation Loop

for t = 0:timeStep:totalTime
    % Update spitz concentration
    spitzArray = spitzArray + timeStep*(spitzDiff*del2(spitzArray) + secretionArray*spitzSecretion - spitzArray*spitzDecay);
     
    % Perform cell logic
    if mod(t,logicStep) == 0
        [cellSpitz, ~, cellSize] = integratedCellIntensity(cellBoundaries, spitzArray);
        divisionOccurs = find(divProb(cellSpitz, cellSize*xStep)>rand(1,length(cellSpitz)));
        [cellBoundaries, engrailed] = proliferationSimulation(cellBoundaries, divisionOccurs, engrailed);
        
        [cellSpitz, ~, cellSize] = integratedCellIntensity(cellBoundaries, spitzArray);
        deathOccurs = find(deathProb(cellSpitz, cellSize*xStep)>rand(1,length(cellSpitz)));
        [cellBoundaries, engrailed] = deathSimulation(cellBoundaries, deathOccurs, engrailed);
        
        cellBoundaries(1) = 1;
        cellBoundaries(arrayLength) = 1;
    end
    
    % Output plot
    if mod(t,outputTime) == 0
        
        subplot(2,1,1)
        [~, cellPlot, ~] = integratedCellIntensity(cellBoundaries, spitzArray);
        plot(x, spitzArray, x, cellPlot, x, cellBoundaries)
        axis([-inf,inf,0,0.6])
        
        subplot(2,1,2)
        imshow(cellImage(cellBoundaries, engrailed), []);
        drawnow
    end
end
