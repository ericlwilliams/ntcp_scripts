function [new_x,new_y,err_y] = rebin(X,Y,new_bins)
% Rebin X to new_bins, Y rebinned too
new_x = zeros(1,length(new_bins));
new_y = zeros(1,length(new_bins));
n_counts = zeros(1,length(new_bins));

for i=1:length(X),
    cur_x = X(i);
    for j=1:length(new_bins),
        if cur_x < new_bins(j),
            new_x(j)=new_x(j)+X(i);%sum for avg
            new_y(j)=new_y(j)+Y(i,1);
            n_counts(j)=n_counts(j)+1;
            break
        end
    end
end
new_x=new_x./n_counts;
new_y=new_y./n_counts;
err_y = (n_counts.^(-.5));
