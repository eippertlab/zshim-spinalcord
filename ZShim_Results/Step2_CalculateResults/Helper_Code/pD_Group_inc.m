function output = pD_Group_inc(x,y)
% calculate percentage difference increases using the group mean
    output =  (mean(y)/mean(x) - 1) * 100;
end