function print_fig(cur_fig,loc,name,type)
    figure_path = strcat(loc,name);
    if strcmp(type,'png'),
        figure_path = strcat(figure_path,'.png');
    end
    print_command = strcat('-d',type);
    disp(['Printing to file: ',figure_path]);
    print(cur_fig,figure_path,print_command);
end