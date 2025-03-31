clear all;close all;clc
num_ele=20;
num_snapshot=10;
dd = 0.5; %d/lambda
loc_array=ULA_gen(num_ele,dd); %生成ULA

snr_mat=[-20:2:10];
Monte=20000;
th=2;
MSE_CBIHT=zeros(1,length(snr_mat));
MSE_iCBIHT=zeros(1,length(snr_mat));
MSE_MUSIC=zeros(1,length(snr_mat));
MSE_MUSIC_arcsine=zeros(1,length(snr_mat));
MSE_MUSIC_unclipped=zeros(1,length(snr_mat));
t_CBIHT=zeros(1,length(snr_mat));
t_iCBIHT=zeros(1,length(snr_mat));
t_MUSIC=zeros(1,length(snr_mat));
t_MUSIC_arcsine=zeros(1,length(snr_mat));
t_MUSIC_unclipped=zeros(1,length(snr_mat));
iter_CBIHT=zeros(1,length(snr_mat));
iter_iCBIHT=zeros(1,length(snr_mat));
normerr_CBIHT=zeros(1,length(snr_mat));
normerr_iCBIHT=zeros(1,length(snr_mat));
for isnr=1:length(snr_mat)
    suc_iCBIHT(isnr)=Monte;
    suc_CBIHT(isnr)=Monte;
    suc_MUSIC(isnr)=Monte;
    suc_MUSIC_arcsine(isnr)=Monte;
    suc_MUSIC_unclipped(isnr)=Monte;
    for iMonte=1:Monte
        doa=[10.3];
        num_source=length(doa);
        snr=snr_mat(isnr);
        xn=gen_data(doa,loc_array,snr,num_snapshot);
        zn=sign_c(xn);
        [doa_CBIHT,p_CBIHT,grid_CBIHT,t_CBIHT_temp,iter_CBIHT_temp,err_CBIHT_temp]=CBIHT(zn,loc_array,num_source);
        [doa_iCBIHT,p_iCBIHT,grid_iCBIHT,t_iCBIHT_temp,iter_iCBIHT_temp,err_iCBIHT_temp]=iCBIHT(zn,loc_array,num_source);
        [doa_MUSIC,p_MUSIC,grid_MUSIC,t_MUSIC_temp]=MUSIC(zn,loc_array,num_source);
        [doa_MUSIC_arcsine,p_MUSIC_arcsine,grid_MUSIC_arcsine,t_MUSIC_arcsine_temp]=MUSIC_arcsine(zn,loc_array,num_source);
        [doa_MUSIC_unclipped,p_MUSIC_unclipped,grid_MUSIC_unclipped,t_MUSIC_unclipped_temp]=MUSIC(xn,loc_array,num_source);
        err_CBIHT=doa-doa_CBIHT;
        err_iCBIHT=doa-doa_iCBIHT;
        err_MUSIC=doa-doa_MUSIC;
        err_MUSIC_arcsine=doa-doa_MUSIC_arcsine;
        err_MUSIC_unclipped=doa-doa_MUSIC_unclipped;
      
        if max(abs(err_CBIHT))>th
            suc_CBIHT(isnr)=suc_CBIHT(isnr)-1;
        end
        if max(abs(err_iCBIHT))>th
            suc_iCBIHT(isnr)=suc_iCBIHT(isnr)-1;
        end        
        if max(abs(err_MUSIC))>th
            suc_MUSIC(isnr)=suc_MUSIC(isnr)-1;
        end        
        if max(abs(err_MUSIC_arcsine))>th
            suc_MUSIC_arcsine(isnr)=suc_MUSIC_arcsine(isnr)-1;
        end        
        if max(abs(err_MUSIC_unclipped))>th
            suc_MUSIC_unclipped(isnr)=suc_MUSIC_unclipped(isnr)-1;
        end       
        
        
        MSE_CBIHT(isnr)=MSE_CBIHT(isnr)+err_CBIHT*err_CBIHT'/Monte/num_source;
        MSE_iCBIHT(isnr)=MSE_iCBIHT(isnr)+err_iCBIHT*err_iCBIHT'/Monte/num_source;
        MSE_MUSIC(isnr)=MSE_MUSIC(isnr)+err_MUSIC*err_MUSIC'/Monte/num_source;
        MSE_MUSIC_arcsine(isnr)=MSE_MUSIC_arcsine(isnr)+err_MUSIC_arcsine*err_MUSIC_arcsine'/Monte/num_source;
        MSE_MUSIC_unclipped(isnr)=MSE_MUSIC_unclipped(isnr)+err_MUSIC_unclipped*err_MUSIC_unclipped'/Monte/num_source;
        t_CBIHT(isnr)=t_CBIHT(isnr)+t_CBIHT_temp/Monte;
        t_iCBIHT(isnr)=t_iCBIHT(isnr)+t_iCBIHT_temp/Monte;
        t_MUSIC(isnr)=t_MUSIC(isnr)+t_MUSIC_temp/Monte;
        t_MUSIC_arcsine(isnr)=t_MUSIC_arcsine(isnr)+t_MUSIC_arcsine_temp/Monte;
        t_MUSIC_unclipped(isnr)=t_MUSIC_unclipped(isnr)+t_MUSIC_unclipped_temp/Monte;
        iter_CBIHT(isnr)=iter_CBIHT(isnr)+iter_CBIHT_temp/Monte;
        iter_iCBIHT(isnr)=iter_iCBIHT(isnr)+iter_iCBIHT_temp/Monte;
        normerr_CBIHT(isnr)=normerr_CBIHT(isnr)+err_CBIHT_temp/Monte;
        normerr_iCBIHT(isnr)=normerr_iCBIHT(isnr)+err_iCBIHT_temp/Monte;
        disp(['进度：',num2str(((isnr-1)*Monte+iMonte)/Monte/length(snr_mat)*100),'%']);
    end
end
RMSE_CBIHT=sqrt(MSE_CBIHT);
RMSE_iCBIHT=sqrt(MSE_iCBIHT);
RMSE_MUSIC=sqrt(MSE_MUSIC);
RMSE_MUSIC_arcsine=sqrt(MSE_MUSIC_arcsine);
RMSE_MUSIC_unclipped=sqrt(MSE_MUSIC_unclipped);



figure;
semilogy(snr_mat,RMSE_MUSIC_unclipped,'-.')
hold on
semilogy(snr_mat,RMSE_MUSIC_arcsine,'--')
semilogy(snr_mat,RMSE_MUSIC,'-.')
semilogy(snr_mat,RMSE_CBIHT,'--')
semilogy(snr_mat,RMSE_iCBIHT)
grid on
xlabel('信噪比（dB）')
ylabel('均方根误差（°）')
legend(['MUSIC（理想含噪观测）'],['1-bit MUSIC（反正弦近似）'],['1-bit MUSIC（直接近似）'],['CBIHT'],['提出的方法'])

figure;
plot(snr_mat,suc_MUSIC_unclipped/Monte,'-.')
hold on
plot(snr_mat,suc_MUSIC_arcsine/Monte,'--')
plot(snr_mat,suc_MUSIC/Monte,'-.')
plot(snr_mat,suc_CBIHT/Monte,'--')
plot(snr_mat,suc_iCBIHT/Monte)
grid on
xlabel('信噪比（dB）')
ylabel('估计成功率')
legend(['MUSIC（理想含噪观测）'],['1-bit MUSIC（反正弦近似）'],['1-bit MUSIC（直接近似）'],['CBIHT'],['提出的方法'])


% %%
figure;
semilogy(snr_mat,t_MUSIC_unclipped,'-.')
hold on
semilogy(snr_mat,t_MUSIC_arcsine,'--')
semilogy(snr_mat,t_MUSIC,'-.')
semilogy(snr_mat,t_CBIHT,'--')
semilogy(snr_mat,t_iCBIHT)
grid on
xlabel('信噪比（dB）')
ylabel('平均运行时间（s）')
legend(['MUSIC（理想含噪观测）'],['1-bit MUSIC（反正弦近似）'],['1-bit MUSIC（直接近似）'],['CBIHT'],['提出的方法'])


figure;
semilogy(snr_mat,iter_CBIHT,'--','color',[0.49,0.18,0.56])
hold on
semilogy(snr_mat,iter_iCBIHT,'color',[0.47,0.67,0.19])
grid on
xlabel('信噪比（dB）')
ylabel('平均迭代次数')
legend(['CBIHT'],['提出的方法'])




