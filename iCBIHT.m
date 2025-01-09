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
S=S-step_size*A'*(sign_c(A*S)-Y); % �ݶ��½�
S_temp=S;
for i=1:K
    S_norm(i,:)=norm(S(i,:)); % �ռ���
end
power=S_norm;
% plot(grid,S_norm)
% hold on
[p,index]=findpeaks(S_norm); % �ҳ���ֵ
[p,index_temp]=sort(p);
S=zeros(K,N);
index_peak=index(index_temp(end-D+1:end)); % ��ֵ���������
S(index_peak,:)=S_temp(index_peak,:); % Ӳ��ֵ
RX=sign_c(A*S);
RX_real=real(RX);RX_imag=imag(RX);
err_real=RX_real.*real(Y);err_imag=RX_imag.*imag(Y);
err_real(err_real>0)=0;err_real=abs(err_real);
err_imag(err_imag>0)=0;err_imag=abs(err_imag);
err=sum(sum(err_real))+sum(sum(err_imag));
err=err/M/N/2;
if iter>=max_iter || err<max_err||min_grid_new<=min_grid % ����������/��С���/��С������
    break;
else % ���ݲ����»�����
    S_peak=S_temp(index_peak,:); % ��ֵ�����ź�ֵ
    S_peak_l=S_temp(index_peak-1,:); % ��ֵ������ź�ֵ
    S_peak_r=S_temp(index_peak+1,:); % ��ֵ�Ҹ����ź�ֵ
    grid_peak=grid(index_peak);
    grid_add_l=grid(index_peak)/2+grid(index_peak-1)/2; % ����Ӹ��ĽǶ�ֵ
    grid_add_r=grid(index_peak)/2+grid(index_peak+1)/2; % ����Ӹ��ĽǶ�ֵ
    [~,index_del]=sort(S_norm); % ɾ�����������
    index_del=index_del(1:2*D);
    del_grid=grid(index_del);
    grid(index_del)=[]; % ɾ���׷���̫С�Ļ�
    grid=sort([grid,grid_add_l,grid_add_r]); % ���º�ĽǶȷ�Χ
    A = exp(-1j*2*pi*loc_array'*sin(grid*pi/180)); % steering matrixv
    for i=1:D
        index_peak(i)=find(grid==grid_peak(i)); % ���º�Ļ������ֵ����
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
t=toc; % ����ʱ��
% p=10*log10(S_norm/max(S_norm));
[peak,index]=findpeaks(p);
[~,index2]=sort(peak);
doa=grid(sort(index(index2(end-D+1:end))));
grid=grid(2:end-1);
p=p(2:end-1);
p=p/max(p);
% p=10*log10(p);
    
    
