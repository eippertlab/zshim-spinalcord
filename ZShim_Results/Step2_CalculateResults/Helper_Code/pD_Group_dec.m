function output = pD_Group_dec(x,y)
% calculate percentage difference decreases using the group mean
    output =  (1- mean(y)/mean(x)) * 100;
end