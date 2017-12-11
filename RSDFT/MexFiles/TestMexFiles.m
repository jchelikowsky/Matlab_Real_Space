function success=TestMexFiles()
% checks for the existence of all the required mex files and tests one to
% determine if all the mex files are compatible with the current version of
% MATLAB

    success=1;
    %check to see that all needed mex files exist
    if (exist('PsuedoDiagLoops','file')~=3)
        success=0;
    elseif (exist('Rayleighritz','file')~=3)
        success=0;
    elseif (exist('malloc','file')~=3)
        success=0;
    elseif (exist('Ceperly_Alder','file')~=3)
        success=0;
    elseif (exist('pcgMexFile','file')~=3)
        success=0;
    elseif (exist('lanWhileLoop','file')~=3)
        success=0;
    elseif (exist('lanWhileLoopForChsubsp','file')~=3)
        success=0;
    elseif (exist('findFirstColumnWithNonZeroElement','file')~=3)
        success=0;
    elseif (exist('sparseMatrixVectorMultiply','file')~=3)
        success=0;
    elseif (exist('BTimesv1MinusbetTimesv0','file')~=3)
        success=0;
    elseif (exist('vMinusalpTimesv1','file')~=3)
        success=0;
    elseif (exist('chebyshevfilterDegree1','file')~=3)
        success=0;
    elseif (exist('chebyshevfilterDegreeN','file')~=3)
        success=0;
    end

    % check for mex file version compatability with MATLAB version.  Assume all mex files were
    % compiled under the same MATLAB version
    try
        tempVar=malloc(1,1);
        clear tempVar;
    catch
        success=0;
    end    

return;

