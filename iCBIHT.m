function [doa,p,grid,t,iter,err]=iCBIHT(xn,loc_array,num_source)
grid=[-90:2:90];

max_iter=100;
min_grid=0.1;
min_grid_new=2;
Y=xn;
D=num_source;
A = exp(-1j*2*pi*loc_array'*sin(grid*pi/180)); % steering matrix
[M,N]=size(xn);
step_size=0.001;
max_err=0.05;
K=length(grid);

%% initialization
S=pinv(A)*Y;
iter=0;
err=inf;

%% iteration
tic
while(1)
iter=iter+1;
S=S-step_size*A'*(sign_c(A*S)-Y); % 梯度下降
S_temp=S;
for i=1:K
    S_norm(i,:)=norm(S(i,:)); % 空间谱
end
power=S_norm;
% plot(grid,S_norm)
% hold on
[p,index]=findpeaks(S_norm); % 找出峰值
[p,index_temp]=sort(p);
S=zeros(K,N);
index_peak=index(index_temp(end-D+1:end)); % 峰值的网格序号
S(index_peak,:)=S_temp(index_peak,:); % 硬阈值
RX=sign_c(A*S);
RX_real=real(RX);RX_imag=imag(RX);
err_real=RX_real.*real(Y);err_imag=RX_imag.*imag(Y);
err_real(err_real>0)=0;err_real=abs(err_real);
err_imag(err_imag>0)=0;err_imag=abs(err_imag);
err=sum(sum(err_real))+sum(sum(err_imag));
err=err/M/N/2;
if iter>=max_iter || err<max_err||min_grid_new<=min_grid % 最大迭代次数/最小误差/最小网格间距
    break;
else % 回溯并更新基向量
    S_peak=S_temp(index_peak,:); % 峰值处的信号值
    S_peak_l=S_temp(index_peak-1,:); % 峰值左格点的信号值
    S_peak_r=S_temp(index_peak+1,:); % 峰值右格点的信号值
    grid_peak=grid(index_peak);
    grid_add_l=grid(index_peak)/2+grid(index_peak-1)/2; % 左添加格点的角度值
    grid_add_r=grid(index_peak)/2+grid(index_peak+1)/2; % 右添加格点的角度值
    [~,index_del]=sort(S_norm); % 删除的网格序号
    index_del=index_del(1:2*D);
    del_grid=grid(index_del);
    grid(index_del)=[]; % 删除谱幅度太小的基
    grid=sort([grid,grid_add_l,grid_add_r]); % 更新后的角度范围
    A = exp(-1j*2*pi*loc_array'*sin(grid*pi/180)); % steering matrixv
    for i=1:D
        index_peak(i)=find(grid==grid_peak(i)); % 更新后的基，其峰值索引
    end
    S=zeros(K,N);
    S(index_peak,:)=S_peak;

%     S(index_peak-2,:)=S_peak_l;
    S(index_peak-1,:)=pinv(A(:,index_peak-1))*A(:,index_peak-2)*S_peak_l;
%     S(index_peak+2,:)=S_peak_r;
    S(index_peak+1,:)=pinv(A(:,index_peak+1))*A(:,index_peak+2)*S_peak_r;
%     for n=1:N
%         S(:,n)=S(:,n)/norm(S(:,n));
%     end
    min_grid_new=min(abs(grid(2:end)-grid(1:end-1)));
%     plot(grid,abs(S))
end
end
for i=1:K
    S_norm(i,:)=norm(S(i,:));
end
p=S_norm;
p=power;
p=[0;p;0];
grid=[-inf,grid,inf];
t=toc; % 运行时间
% p=10*log10(S_norm/max(S_norm));
[peak,index]=findpeaks(p);
[~,index2]=sort(peak);
doa=grid(sort(index(index2(end-D+1:end))));
grid=grid(2:end-1);
p=p(2:end-1);
p=p/max(p);
% p=10*log10(p);
    
    
