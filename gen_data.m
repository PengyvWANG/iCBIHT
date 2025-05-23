function xn=gen_data(doa,loc_array,snr,num_snapshot)
num_ele = length(loc_array); % 阵元数
num_source = length(doa); % 信源个数
A = exp(-1j*2*pi*loc_array'*sin(doa*pi/180)); % steering matrix
source=(randn(num_source,num_snapshot)+1j*randn(num_source,num_snapshot))/sqrt(2)*sqrt(10^(snr/10));
x=A*source;
noise=(randn(num_ele,num_snapshot)+1j*randn(num_ele,num_snapshot))/sqrt(2);
% xn=awgn(x,snr,'measured');
xn=x+noise;