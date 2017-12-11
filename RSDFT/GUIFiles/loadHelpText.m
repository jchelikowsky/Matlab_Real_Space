function cellArrayOfHelp = loadHelpText(helpTopic)

%remove spaces
helpTopic(helpTopic==' ')=[];

path(path,'HelpFiles');

filename=horzcat(helpTopic,'.txt');

cellArrayOfHelp={};

try
    cellArrayOfHelp=textread(filename,'%s','delimiter', '\n','whitespace', '');
catch
    exceptionID=lasterror();

    if (findstr(exceptionID.message,'File not found'))
        warndlg({'Help File Not Found',horzcat('File name is: ',filename)},'Exception Thrown');       
    else
        warndlg({'An error occurred while trying to read help file',horzcat('ERROR MESSAGE: ',exceptionID.message)},'Exception Thrown');  
    end    
    
end    
    
return;