function loc_array=ULA_gen(num_ele,dd)
if rem(num_ele,2)==0
    loc_array = [-num_ele/2*dd+dd/2:dd:num_ele/2*dd-dd/2]; % ��һ����Ԫλ��
else
    loc_array = [-(num_ele-1)/2*dd:dd:(num_ele-1)/2*dd]; % ��һ����Ԫλ��
end