function [ outputArray ] = periodicBoundaryCondition( inputArray )

    outputArray = [inputArray(end) inputArray inputArray(1)];

end

