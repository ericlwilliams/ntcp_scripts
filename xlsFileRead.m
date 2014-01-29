function xlsData=xlsFileRead(xlsFile)
% read all the sheets in the file
    % check if it is a .xls compatable file
    % xlsinfo1 - 'Microsoft Excel Spreadsheet' xlsinfo2 - names of
    % worksheets
    [xlsinfo1,xlsinfo2]=xlsfinfo(xlsFile); 
    
    if ~strcmp('Microsoft Excel Spreadsheet',xlsinfo1)
        xlsData=[]; disp('not a .xls file'); return;
    end
    
    % read each sheet to a structure
    num=length(xlsinfo2);
    xlsData=repmat(struct('SheetName','','xlsRaw',''),[num,1]);
    % xlsData - array of structures with properties:
    %   sheetName: 'Whole', 'Ipsi', 'Contra', etc...
    %   xlsRaw: raw data from sheet (as returned by xlsread below)
    for k=1:num
        [~,~,raw]=xlsread(xlsFile,xlsinfo2{k});
        xlsData(k).SheetName=xlsinfo2{k};
        xlsData(k).xlsRaw=raw;
    end
