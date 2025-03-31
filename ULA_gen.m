function loc_array=ULA_gen(num_ele,dd)
if rem(num_ele,2)==0
    loc_array = [-num_ele/2*dd+dd/2:dd:num_ele/2*dd-dd/2]; % 归一化阵元位置
else
    loc_array = [-(num_ele-1)/2*dd:dd:(num_ele-1)/2*dd]; % 归一化阵元位置
end