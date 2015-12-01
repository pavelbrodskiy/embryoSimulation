% This class represents a 1D sheet of cells in order of adjacency

classdef SimulationCellArray
    properties
        cellArray@SimulationCell       % Array containing the cell sheet
    end
    methods
        % Contructor
        function obj = SimulationCellArray(ACells, PCells, xStep)
            % Initialize anterior compartment (left side)
            obj.cellArray(1:ACells) = SimulationCell(false, xStep);
            
            % Initialize posterior compartment (right side)
            obj.cellArray((ACells+1):(ACells+PCells)) = SimulationCell(true, xStep);
        end
        
        % Return pixel volumes
        function volumes = getVolumes(obj)
            for cell = length(obj.cellArray):-1:1
                volumes(cell) = length(obj.cellArray(cell).spitzArray);
            end
        end
        
        % Return ligand arrays
        function [spitzArray, winglessArray] = ligandArrays(obj)
            spitzArray = [obj.cellArray(:).spitzArray];
            winglessArray = [obj.cellArray(:).winglessArray];
        end
        
        % Divide array: for each cell where the division timer has passed
        % the simulation time, duplicate the cell and set both daughter
        % cells to half the volume of the mother cell. Divide the
        % underlying ligand concentrations into two equal domains.
        function obj = divideCellArray(obj, t)
            for cell = length(obj.cellArray):-1:1
                timer = obj.cellArray(cell).divisionTimer;
                if ~isempty(timer) && timer(1) <= t
                    % Remove the division from the programming
                    timer(1) = [];
                    
                    % Create the daughter cells
                    cellTemp1 = obj.cellArray(cell);
                    cellTemp1.divisionTimer = timer;
                    cellTemp1.volume = cellTemp1.volume/2;
                    cellTemp2 = cellTemp1;
                    
                    firstIndex = uint16(length(cellTemp1.spitzArray)/2+.5);
                    cellTemp2.spitzArray = cellTemp1.spitzArray(firstIndex:end);
                    cellTemp2.winglessArray = cellTemp1.winglessArray(firstIndex:end);
                    cellTemp1.spitzArray = cellTemp1.spitzArray(1:(firstIndex-1));
                    cellTemp1.winglessArray = cellTemp1.winglessArray(1:(firstIndex-1));
                    
                    % Replace the mother cell with the daughter cells
                    obj.cellArray = [obj.cellArray(1:(cell-1)) cellTemp1 cellTemp2 obj.cellArray((cell+1):end)];
                end
            end
        end
        
        % Grow cell array: Determine which cells grow at this timepoint
        % as a function of mechanical and chemical properties, then perform
        % growth and apoptosis.
        function obj = growKillCellArray(obj, deathProb, growthRate, deathSpeed, dt)
            for cell = length(obj.cellArray):-1:1
                timer = obj.cellArray(cell).divisionTimer;
                
                spitz = mean(obj.cellArray(cell).spitzArray);
                wg = mean(obj.cellArray(cell).winglessArray);
                size = obj.cellArray(cell).volume;
                en = obj.cellArray(cell).engrailed;
                
                %if isempty(timer)
                    if deathProb(spitz,wg,en,size)*dt > rand(1)
                        obj.cellArray(cell).dying = true;
                    end
                %end
                
                if obj.cellArray(cell).dying
                    growth = -deathSpeed*dt;
                else
                    growth = growthRate(spitz,wg,en,size)*dt;
                end
                obj.cellArray(cell).volume = obj.cellArray(cell).volume + growth;
                newDomainSize = obj.cellArray(cell).volume/obj.cellArray(cell).xStep;
                if newDomainSize <= 0
                    obj.cellArray(cell) = [];
                else
                    initialSpi = sum(obj.cellArray(cell).spitzArray);
                    initialWg = sum(obj.cellArray(cell).winglessArray);
                    
                    obj.cellArray(cell).spitzArray = imresize(obj.cellArray(cell).spitzArray,[1,newDomainSize]);
                    obj.cellArray(cell).winglessArray = imresize(obj.cellArray(cell).winglessArray,[1,newDomainSize]);
                    
                    finalSpi = sum(obj.cellArray(cell).spitzArray);
                    finalWg = sum(obj.cellArray(cell).winglessArray);
                    
                    obj.cellArray(cell).spitzArray = obj.cellArray(cell).spitzArray * initialSpi / finalSpi;
                    obj.cellArray(cell).winglessArray = obj.cellArray(cell).winglessArray * initialWg / finalWg;
                end
                
            end
        end
        
        
        % Find secretion function: This returns a function for how ligands
        % are secreted (where 1 represents maximum secretion of ligand).
        % Currently this is a step function with 1 on secreting cells and 0
        % on non-secreting cells.
        function [spitzSecretionArray, winglessSecretionArray, PCompartment] = findSecretionSites(obj, spitzCells, winglessCells)
            
            [spitzArray, ~] = ligandArrays(obj);
            domainSize = length(spitzArray);
            
            spitzSecretionArray = zeros(1, domainSize);
            winglessSecretionArray = spitzSecretionArray;
            PCompartment = spitzSecretionArray;
            firstPCell = [];
            
            for cell = 1:length(obj.cellArray)
                if obj.cellArray(cell).engrailed == true && isempty(firstPCell)
                    firstPCell = cell;
                end
            end
            
            xStep = obj.cellArray(1).xStep;
            
            % Find P compartment nodes
            startOfP = 1;
            if firstPCell > 1
                for cell = 1:(firstPCell-1)
                    startOfP = startOfP + length(obj.cellArray(cell).spitzArray);
                end
                
                
                PCompartment(startOfP:end) = true;
                
                % Find spitz secreting cells (left side of A compartment)
                divisionNumber = length(spitzCells) - length(obj.cellArray(1).divisionTimer);
                endOfSpitz = 0;
                for cell = 1:spitzCells(divisionNumber)
                    if ~obj.cellArray(cell).engrailed
                        endOfSpitz = endOfSpitz + length(obj.cellArray(cell).spitzArray);
                    end
                end
                
                if endOfSpitz > 0
                    spitzSecretionArray(1:endOfSpitz) = 1;
                end
                
                % Find wingless secreting cells (right side of A compartment)
                
                divisionNumber = length(spitzCells) - length(obj.cellArray(firstPCell-1).divisionTimer);
                firstWinglessCell = max([1, (firstPCell-winglessCells(divisionNumber))]);
                
                startofWg = 1;
                if firstWinglessCell > 1
                    for cell = 1:(firstWinglessCell-1)
                        startofWg = startofWg + length(obj.cellArray(cell).spitzArray);
                    end
                end
                
                endofwg = startofWg - 1;
                for cell = firstWinglessCell:(firstPCell-1)
                    endofwg = endofwg + length(obj.cellArray(cell).spitzArray);
                end
                
                winglessSecretionArray(startofWg:endofwg) = 1;
                
            end
        end
        
        % Update concentrations: This resets the concentration array in
        % each cell to synchronize with the concentrations in the domain.
        function obj = updateConcentrations(obj, spitzArray, winglessArray)
            firstIndex = 1;
            for cell = 1:length(obj.cellArray)
                cellSize = length(obj.cellArray(cell).spitzArray);
                obj.cellArray(cell).spitzArray = spitzArray(firstIndex:(firstIndex+cellSize-1));
                obj.cellArray(cell).winglessArray = winglessArray(firstIndex:(firstIndex+cellSize-1));
                firstIndex = firstIndex + cellSize;
            end
        end
        
        % Get dying cells: Returns a logical array containing the
        % identities of cells undergoing apoptosis.
        function dyingCells = getDyingCells(obj)
            dyingCells = [obj.cellArray(:).dying];
        end
    end
end