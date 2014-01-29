function [dDVH,cDVH,ddvh_var_name,cdvh_var_name]=import_ddvh_file(fileToRead1,fileNum)
%IMPORTFILE(FILETOREAD1)
%  Imports data from the specified file
%  FILETOREAD1:  file to read

%  Auto-generated by MATLAB on 01-Aug-2012 10:25:51

% Import the file
newData1 = importdata(fileToRead1);
clear dDVH;
% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    var_name = char(vars{i});
    %assignin('base', vars{i}, newData1.(vars{i}));
    if strfind(char(vars{i}),'textdata')
        continue;
    end
   %assignin('base', cell(strcat(char(vars{i}),char(fileNum))), newData1.(vars{i}));
   if fileNum<10,
      tmp_num = int2str(fileNum);
      fileNumStr=strcat('0',tmp_num);
   else
       fileNumStr=int2str(fileNum);      
   end
   ddvh_var_name = strcat('dDVH',fileNumStr);
   dDVH=newData1.(vars{i});
   assignin('base', ddvh_var_name, dDVH);
   
   % make cDVH
   cdvh_var_name = strcat('cDVH',fileNumStr);
   cDVH = cDVH_from_dDVH(dDVH);
   assignin('base',cdvh_var_name,cDVH);
   end
end
