function [] = logError(exception)
% takes in struct from lasterror() and saves important info to errorlog.txt
globalVariables=whos('global');

variableFound=0;
for i=1:size(globalVariables,1)
    if (strcmp(globalVariables(i).name,'enableErrorLogging'))
        variableFound=1;
    end
end


if (variableFound==1)
    eval(horzcat('global ','enableErrorLogging'));
else 
    enableErrorLogging=0;
end

if (variableFound==0 || enableErrorLogging==1)
    report=exception.message;
    stack=exception.stack;


    fileID=fopen('errorlog.txt','a');


    fprintf(fileID,'//////////////////////////////////////////////////////////////////\n');

    dateAndTime=clock();
    dateString=datestr(dateAndTime,21);


    fprintf(fileID,'%s\n\n',dateString);

    [compStr, maxSize, endian]=computer;
    fprintf(fileID,'Computer type: %s    wordSize: %d    endian: %s \n\n', compStr,maxSize,endian);

    fprintf(fileID,'MATLAB version= %s\n\n',version());

    fprintf(fileID,'%s\n\n',report);

    
    for i=1:max(size(stack))
        fprintf(fileID,'file:%s, line #:%d\n\n',stack(i).file,stack(i).line);
    end



    fprintf(fileID,'Global Variables\n');
    variableValue=-1;
    for i=1:size(globalVariables,1)
        eval(horzcat('global ',globalVariables(i).name));
        eval(horzcat('variableValue=',globalVariables(i).name,';'));
        fprintf(fileID,'%s=%f\n',globalVariables(i).name,variableValue);
    end


    fprintf(fileID,'\nSaved rsdft.out\n\n');

    try
        cellArrayOfrsdft=textread('rsdft.out','%s','delimiter', '\n','whitespace', '');
    catch
        exceptionID=lasterror();

        if (findstr(exceptionID.message,'File not found'))
            warndlg({'Help File Not Found',horzcat('File name is: ',filename)},'Exception Thrown');
        else
            warndlg({'An error occurred while trying to read help file',horzcat('ERROR MESSAGE: ',exceptionID.message)},'Exception Thrown');
        end

    end

    for i=1:size(cellArrayOfrsdft,1)
        fprintf(fileID,'        %s\n',cellArrayOfrsdft{i});
    end


    fprintf(fileID,'//////////////////////////////////////////////////////////////////\n');

    fclose(fileID);
end
