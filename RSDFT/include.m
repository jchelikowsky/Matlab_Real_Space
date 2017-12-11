%%
% Basic Settings

global CG_prec poldeg diagmeth adaptiveScheme enableChargeDensityVisualization
global OPTIMIZATIONLEVEL;

% inlcude paths here so that main.m can also run on its own!
path(path,'MIXER');
path(path,'SPLINES');
path(path,'GUIFiles');
path('DiagonalizationFiles',path); %added to top of path because MATLAB also has a function named pcg
path(path,'IonicPotentialFindingFiles');
path(path,'MexFiles');
path(path,'LaplacianMatrixFiles');

%-------------------- whether or not to precondition CG
%
CG_prec = 0;
%
%-------------------- polynomial degree for chebyshev
%
poldeg = 10;
%
%-------------------- method for diagonalization
%  diagmeth == 0 --> Lanczos 1st step and chebyshev filtering thereafter
%  diagmeth == 1 --> Lanczos all the time
%  diagmeth == 2 --> Full-Chebyshev subspace iteration first step
%                     chebyshev filtering thereafter.
%  diagmeth == 3 --> First step filtering of random initial vectors,
%                   (this uses less memory and faster than diagonalization)
%                     then chebyshev subspace filtering thereafter
%diagmeth = 0;
diagmeth = 3;
% adaptiveScheme == 0    Do not use an adaptive scheme
% adaptiveScheme == 1    Allow for the changing of parameters used by
%                         lanczos and chefsi1 to increase speed,  slight risk of longer time to
%                         converge
adaptiveScheme=0;
%
% control whether or not the visualization of the charge density is
% enabled or disabled
% 0 is disabled
% 1 is enabled
enableChargeDensityVisualization=0;

if (isunix()==1)
    enableChargeDensityVisualization=0;
end    
%
% control which user interface is used
% 0 is text based
% not 0 is graphical user interface
userInterfaceControl=1;

% 0= Only use MATLAB code
% 1= Use some precompiled c code
OPTIMIZATIONLEVEL=0;

%%
%%-------------------------------
% Advanced Settings
global fd_order maxits tol Fermi_temp
%%%
%%%%Defaut Technical Parameters
%%%%
fd_order = 8;                         %% order of finite difference scheme {8}
maxits   = 40;                         %% max SCF iterations               {40}
tol      = 1.e-03;                     %% tolerance for SCF iteration.     {1.e-03}
Fermi_temp  =  500.0 ;                %%  Smear out fermi level            {500.0}


%% debugging and advanced profiling settings
global enableMexFilesTest enableErrorLogging profileDataOutput
global enableOSCheckForVisualization
% when set to 1, both mex files and MATLAB code are run and are then
% compared for accuracy, if accuracy is outside bounds, error is logged
% OPTIMIZATIONLEVEL must be set to 0
% turning this option on slows down the code tremendously

% tested files- BTimesv1MinusbetTimesv0, Ceperly_Alder,
% chebyshevfilterDegree1, chebyshevfilterDegreeN, pcgMexFile, Rayleighritz,
% vMinusalpTimesv1, psuedoDiag
enableMexFilesTest=0;

if (enableMexFilesTest==1)
    OPTIMIZATIONLEVEL=0;
end    

% when set to 1, when an error occurs, relevent information is saved to the
% file, errorlog.txt
enableErrorLogging=1;

% Set profileDataOutput to 0 writes to command window
% when set to 1, writes profiling code to file, file name is date and time
% of profile .txt
% any thing else is no output

profileDataOutput=-1;


% if 0, does not check the current os when running the charge density
% visualiztaion, the os is checked because the visualization can cause
% crashes on Linux machines
% The cause of the crashing seems to be the use of opengl for rendering so
% opengl is disabled, if problems do arise on unix machines even without
% the use of opengl, set this value to 1 and the use of the charge density
% visualization on unix machines will be disabled
enableOSCheckForVisualization=0;
