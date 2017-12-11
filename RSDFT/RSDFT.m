% a script that starts RSDFT
% its purpose is to choose the correct
% user interface based on the current MATLAB version

% also adds paths for MATLAB to look for files in
clear all;





% tell MATLAB not to use opengl for rendering
% opengl can cause crashes on UNIX machines
%opengl('neverselect');

include;

versionNumber=version;
%get the first 2 digits of the version number, include decimal point
versionNumber=versionNumber(1:3);

% if the verion of matlab is 7.0.* or earlier, or is explicitly selected,
% use the text based user interface
if (userInterfaceControl==0 || str2double(versionNumber)<=7.0)
    if (OPTIMIZATIONLEVEL~=0 && TestMexFiles()==0)
        disp('There is an error with the mex files, using OPTIMIZATIONLEVEL 0');
        OPTIMIZATIONLEVEL=0;
    end
    try
        main;
    catch
        exceptionID=lasterror();
        
        logError(exceptionID);
        
        disp('An error occurred during execution.')
        disp(horzcat('ERROR MESSAGE: ',exceptionID.message))
        disp('Please check your settings and if this problem persists, please contact the proper people.');
    end    
else
    MoleculeCreationGUI;
end