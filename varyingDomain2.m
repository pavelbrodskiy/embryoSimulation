%% Definitions
% Simulation Parameters
domainSize = 200;               %microns
xStep = 0.1;                    %microns
totalTime = 1000;               %s
timeStep = 0.001;               %s
logicStep = 0.1;                %s
outputTime = 0.1;               %s

% Spitz Parameters
spitzSecretion = 1;             %uM
spitzDecayA = 0.05;             %Percent decayed per second
spitzDecayP = 0.05;             %Percent decayed per second
spitzDiff = 222;                %uM/s^2

% Cell Parameters
cellSize = 5;                   %microns
deathProb = @(spitz,size) 0;
timeToDivision = 5;             %s

%% Initialization
% Initialize Cells
Cell.size = cellSize;
Cell.divisionsLeft = 2;
Cell.divisionTimer = timeToDivision;
Cell.engrailed = 1;

PCells = mod(1:16,3);
cellArray(PCells) = Cell;

Cell.engrailed = 0;
cellArray(~PCells) = Cell;

%Initialize Domain
x = 0:xStep:domainSize;
arrayLength = length(x);
secretionArray = zeros(1,arrayLength);
spitzArray = secretionArray;



%% Simulation Loop
for t = 0:timeStep:totalTime
    %Generate Domain Maps
    clear PDomain ADomain 
    tempPos = 100/xStep;
    engrailed = [cellArray.engrailed];
    for i = 1:length(cellArray)
        if engrailed(i)
            PDomain(tempPos:(tempPos + cellArray.size)) = 1;
            ADomain(tempPos:(tempPos + cellArray.size)) = 0;
        else
            ADomain(tempPos:(tempPos + cellArray.size)) = 1;
            PDomain(tempPos:(tempPos + cellArray.size)) = 0;
        end
        tempPos = tempPos + cellArray.size;
    end
    
    spitzDecay = PDomain .* spitzDecayP + ADomain .* spitzDecayA;
    
    %Update Ligand Concentrations
    spitzArray = spitzArray + timeStep * (spitzDiff * del2(spitzArray) - spitzDecay .* spitzArray);
    
    
    
    
    
    
    
    
end








