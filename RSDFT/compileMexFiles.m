function []=compileMexFiles(enableProgressbar)


if (nargin==0)
    enableProgressbar=0;
end
    

if (enableProgressbar==1)
    h=waitbar(0,'Compiling Mex Files');
end

[arch integerSize]=computer;

if (integerSize<=2^32)
    % code for compiling .mex files

    mex('COPTIMFLAGS=-O3 -w','-outdir','IonicPotentialFindingFiles', horzcat('IonicPotentialFindingFiles',filesep,'PsuedoDiagLoops.c'));
    if (enableProgressbar==1)
        waitbar(.2,h);
    end


    mex('COPTIMFLAGS=-O3 -w','-outdir', 'DiagonalizationFiles',horzcat('DiagonalizationFiles',filesep,'Rayleighritz.c'));
    mex('COPTIMFLAGS=-O3 -w','-outdir', 'MexFiles',horzcat('MexFiles',filesep,'malloc.c'));
    if (enableProgressbar==1)
        waitbar(.4,h);
    end

    mex('COPTIMFLAGS=-O3 -w', '-outdir', 'DiagonalizationFiles',horzcat('DiagonalizationFiles',filesep,'pcgMexFile.c'));
    mex('COPTIMFLAGS=-O3 -w', '-outdir', 'DiagonalizationFiles',horzcat('DiagonalizationFiles',filesep,'Ceperly_Alder.c'));
    if (enableProgressbar==1)
        waitbar(.5,h);
    end

    mex('COPTIMFLAGS=-O3 -w','-outdir', 'DiagonalizationFiles',horzcat('DiagonalizationFiles',filesep,'lanWhileLoop.c'));
    mex('COPTIMFLAGS=-O3 -w','-outdir', 'DiagonalizationFiles',horzcat('DiagonalizationFiles',filesep,'lanWhileLoopForChsubsp.c'));
    if (enableProgressbar==1)
        waitbar(.6,h);
    end

    mex('COPTIMFLAGS=-O3 -w', '-outdir', 'MexFiles',horzcat('MexFiles',filesep,'findFirstColumnWithNonZeroElement.c'));
    mex('COPTIMFLAGS=-O3 -w', '-outdir', 'MexFiles',horzcat('MexFiles',filesep,'sparseMatrixVectorMultiply.c'));
    if (enableProgressbar==1)
        waitbar(.7,h);
    end

    mex('COPTIMFLAGS=-O3 -w','-outdir','DiagonalizationFiles', horzcat('DiagonalizationFiles',filesep,'BTimesv1MinusbetTimesv0.c'));
    if (enableProgressbar==1)
        waitbar(.8,h);
    end

    mex('COPTIMFLAGS=-O3 -w', '-outdir','DiagonalizationFiles', horzcat('DiagonalizationFiles',filesep,'vMinusalpTimesv1.c'));
    if (enableProgressbar==1)
        waitbar(.9,h);
    end

    mex('COPTIMFLAGS=-O3 -w', '-outdir','DiagonalizationFiles',horzcat('DiagonalizationFiles',filesep,'chebyshevfilterDegree1.c'));
    mex('COPTIMFLAGS=-O3 -w', '-outdir','DiagonalizationFiles',horzcat('DiagonalizationFiles',filesep,'chebyshevfilterDegreeN.c'));
    if (enableProgressbar==1)
        waitbar(1,h,'Completed Compilation');
        close(h);
    end
else
   % code for compiling .mex files

    mex('COPTIMFLAGS=-O3 -w','-largeArrayDims','-outdir','IonicPotentialFindingFiles', horzcat('IonicPotentialFindingFiles',filesep,'PsuedoDiagLoops.c'));
    if (enableProgressbar==1)
        waitbar(.2,h);
    end


    mex('COPTIMFLAGS=-O3 -w','-largeArrayDims','-outdir', 'DiagonalizationFiles',horzcat('DiagonalizationFiles',filesep,'Rayleighritz.c'));
    mex('COPTIMFLAGS=-O3 -w','-largeArrayDims','-outdir', 'MexFiles',horzcat('MexFiles',filesep,'malloc.c'));
    if (enableProgressbar==1)
        waitbar(.4,h);
    end

    mex('COPTIMFLAGS=-O3 -w','-largeArrayDims', '-outdir', 'DiagonalizationFiles',horzcat('DiagonalizationFiles',filesep,'pcgMexFile.c'));
    mex('COPTIMFLAGS=-O3 -w','-largeArrayDims', '-outdir', 'DiagonalizationFiles',horzcat('DiagonalizationFiles',filesep,'Ceperly_Alder.c'));
    if (enableProgressbar==1)
        waitbar(.5,h);
    end

    mex('COPTIMFLAGS=-O3 -w','-largeArrayDims','-outdir', 'DiagonalizationFiles',horzcat('DiagonalizationFiles',filesep,'lanWhileLoop.c'));
    mex('COPTIMFLAGS=-O3 -w','-largeArrayDims','-outdir', 'DiagonalizationFiles',horzcat('DiagonalizationFiles',filesep,'lanWhileLoopForChsubsp.c'));
    if (enableProgressbar==1)
        waitbar(.6,h);
    end

    mex('COPTIMFLAGS=-O3 -w','-largeArrayDims', '-outdir', 'MexFiles',horzcat('MexFiles',filesep,'findFirstColumnWithNonZeroElement.c'));
    mex('COPTIMFLAGS=-O3 -w','-largeArrayDims', '-outdir', 'MexFiles',horzcat('MexFiles',filesep,'sparseMatrixVectorMultiply.c'));
    if (enableProgressbar==1)
        waitbar(.7,h);
    end

    mex('COPTIMFLAGS=-O3 -w','-largeArrayDims','-outdir','DiagonalizationFiles', horzcat('DiagonalizationFiles',filesep,'BTimesv1MinusbetTimesv0.c'));
    if (enableProgressbar==1)
        waitbar(.8,h);
    end

    mex('COPTIMFLAGS=-O3 -w','-largeArrayDims', '-outdir','DiagonalizationFiles', horzcat('DiagonalizationFiles',filesep,'vMinusalpTimesv1.c'));
    if (enableProgressbar==1)
        waitbar(.9,h);
    end

    mex('COPTIMFLAGS=-O3 -w','-largeArrayDims', '-outdir','DiagonalizationFiles',horzcat('DiagonalizationFiles',filesep,'chebyshevfilterDegree1.c'));
    mex('COPTIMFLAGS=-O3 -w','-largeArrayDims', '-outdir','DiagonalizationFiles',horzcat('DiagonalizationFiles',filesep,'chebyshevfilterDegreeN.c'));
    if (enableProgressbar==1)
        waitbar(1,h,'Completed Compilation');
        close(h);
    end 
end