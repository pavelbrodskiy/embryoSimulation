% This class represents an individual cell in the sheet

classdef SimulationCell
    properties
        volume@double              % Size of the cell [um]
        spitzArray@double          % Concentration of spitz [AU]
        winglessArray@double       % Concentration of wingless [AU]
        divisionTimer@double       % Time until each division [s]
        engrailed@logical          % Is engrailed turned on? [boolean]
        dying@logical              % Is cell undergoing apoptosis? [boolean]
        xStep@double               % Spatial dimensions of domain [um]
    end
    methods
        % Default constructor: Minimum inputs to create a cell
        function obj = SimulationCell(engrailed, xStep)
            % Values for default constructor
            defaultVolume = 10;
            defaultDivisionTimes = [10000];
            
            % Construct cell from default values
            obj.volume = defaultVolume;
            obj.spitzArray = ones(1,uint16(obj.volume/xStep))*1;
            obj.winglessArray = obj.spitzArray;
            obj.divisionTimer = defaultDivisionTimes;
            obj.engrailed = engrailed;
            obj.dying = false;
            obj.xStep = xStep;
        end
    end
end