load('inoue_transmission_curves.mat')

load('z4_analysis.mat')

temp = log10(t_grid(:,1:end-1)');
temp(temp == -Inf) = 0;
temp2 = trvals(1:end-1);
temp2 = temp2/max(temp2);

mesh(wavs,temp2,temp);hold on
colormap(magma(1000))

plot3(lam_obs,transmission_curves_allz(:,21),ones(size(lam_obs))*max(max(temp))*1.2,'w-','linewidth',2)
axis([3800,9000,0,1,min(min(temp)),max(max(temp))*1.5])
xlabel('\lambda [A]','fontsize',18)
ylabel('T^{IGM}_{\lambda ,z}','fontsize',18)
title('IGM transmission at z=4 [White line: Inoue(2014)]','fontsize',18)


%--------------------------------------------------------------------------

for i = 1:numel(wavs)
    temp_tarr = temp(:,i);
    temp_tarr_smoothed = smooth(temp_tarr,0.1,'loess');
    temp_smoothed(:,i) = temp_tarr_smoothed;
    [maxval(i),maxin(i)] = min(abs(cumsum(10.^temp_tarr)/max(cumsum(10.^temp_tarr))-0.5));
    [maxval2(i),maxin2(i)] = max(temp_tarr_smoothed);
end

%plot(temp(:,100));hold on
%plot(smooth(temp(:,100),0.1,'loess'))



%mesh(wavs,temp2,temp_smoothed);hold on
mesh(wavs,temp2,temp);hold on
colormap(magma(1000))
plot3(wavs,temp2(maxin),ones(size(wavs))*max(max(temp))*1.1,'-','linewidth',2)
plot3(wavs,smooth(temp2(maxin2),0.01,'loess'),ones(size(wavs))*max(max(temp))*1.15,'-','linewidth',2)
plot3(wavs,smooth(temp2(maxin2),0.01,'loess'),ones(size(wavs))*max(max(temp))*1.15,'-','linewidth',2)
%plot3(wavs,smooth(temp2(maxin2),0.1,'rloess'),ones(size(wavs))*max(max(temp))*1.18,'-','linewidth',2)
plot3(lam_obs,transmission_curves_allz(:,21),ones(size(lam_obs))*max(max(temp))*1.2,'w-','linewidth',2)
axis([3800,9000,0,1,min(min(temp)),max(max(temp))*1.5])
xlabel('\lambda [A]','fontsize',18)
ylabel('T^{IGM}_{\lambda ,z}','fontsize',18)
title('IGM transmission at z=4 [White line: Inoue(2014)]','fontsize',18)


%--------------------------------------------------------------------------
% repeat analysis at z = 3.5

load('z3.5_analysis.mat')


temp = log10(t_grid(:,1:end-1)');
temp(temp == -Inf) = 0;
temp2 = trvals(1:end-1);
temp2 = temp2/max(temp2);

mesh(wavs,temp2,temp);hold on
colormap(magma(1000))

plot3(lam_obs,transmission_curves_allz(:,16),ones(size(lam_obs))*max(max(temp))*1.2,'w-','linewidth',2)
axis([3800,9000,0,1,min(min(temp)),max(max(temp))*1.5])
xlabel('\lambda [A]','fontsize',18)
ylabel('T^{IGM}_{\lambda ,z}','fontsize',18)
title('IGM transmission at z=3.5 [White line: Inoue(2014)]','fontsize',18)



%--------------------------------------------------------------------------

for i = 1:numel(wavs)
    temp_tarr = temp(:,i);
    temp_tarr_smoothed = smooth(temp_tarr,0.1,'loess');
    temp_smoothed(:,i) = temp_tarr_smoothed;
    [maxval(i),maxin(i)] = min(abs(cumsum(10.^temp_tarr)/max(cumsum(10.^temp_tarr))-0.5));
    [maxval2(i),maxin2(i)] = max(temp_tarr_smoothed);
end

%plot(temp(:,100));hold on
%plot(smooth(temp(:,100),0.1,'loess'))



%mesh(wavs,temp2,temp_smoothed);hold on
mesh(wavs,temp2,temp);hold on
colormap(magma(1000))
plot3(wavs,temp2(maxin),ones(size(wavs))*max(max(temp))*1.1,'-','linewidth',2)
plot3(wavs,smooth(temp2(maxin2),0.01,'loess'),ones(size(wavs))*max(max(temp))*1.15,'-','linewidth',2)
plot3(wavs,smooth(temp2(maxin2),0.01,'loess'),ones(size(wavs))*max(max(temp))*1.15,'-','linewidth',2)
%plot3(wavs,smooth(temp2(maxin2),0.1,'rloess'),ones(size(wavs))*max(max(temp))*1.18,'-','linewidth',2)
plot3(lam_obs,transmission_curves_allz(:,16),ones(size(lam_obs))*max(max(temp))*1.2,'w-','linewidth',2)
axis([3800,9000,0,1,min(min(temp)),max(max(temp))*1.5])
xlabel('\lambda [A]','fontsize',18)
ylabel('T^{IGM}_{\lambda ,z}','fontsize',18)
title('IGM transmission at z=3.5 [White line: Inoue(2014)]','fontsize',18)




%--------------------------------------------------------------------------
% repeat analysis at z = 4.5

load('z4.5_analysis.mat')


temp = log10(t_grid(:,1:end-1)');
temp(temp == -Inf) = 0;
temp2 = trvals(1:end-1);
temp2 = temp2/max(temp2);

mesh(wavs,temp2,temp);hold on
colormap(magma(1000))

plot3(lam_obs,transmission_curves_allz(:,26),ones(size(lam_obs))*max(max(temp))*1.2,'w-','linewidth',2)
axis([3800,9000,0,1,min(min(temp)),max(max(temp))*1.5])
xlabel('\lambda [A]','fontsize',18)
ylabel('T^{IGM}_{\lambda ,z}','fontsize',18)
title('IGM transmission at z=4.5 [White line: Inoue(2014)]','fontsize',18)

%--------------------------------------------------------------------------

for i = 1:numel(wavs)
    temp_tarr = temp(:,i);
    temp_tarr_smoothed = smooth(temp_tarr,0.1,'loess');
    temp_smoothed(:,i) = temp_tarr_smoothed;
    [maxval(i),maxin(i)] = min(abs(cumsum(10.^temp_tarr)/max(cumsum(10.^temp_tarr))-0.5));
    [maxval2(i),maxin2(i)] = max(temp_tarr_smoothed);
end

%plot(temp(:,100));hold on
%plot(smooth(temp(:,100),0.1,'loess'))



%mesh(wavs,temp2,temp_smoothed);hold on
mesh(wavs,temp2,temp);hold on
colormap(magma(1000))
plot3(wavs,temp2(maxin),ones(size(wavs))*max(max(temp))*1.1,'-','linewidth',2)
plot3(wavs,smooth(temp2(maxin2),0.01,'loess'),ones(size(wavs))*max(max(temp))*1.15,'-','linewidth',2)
plot3(wavs,smooth(temp2(maxin2),0.01,'loess'),ones(size(wavs))*max(max(temp))*1.15,'-','linewidth',2)
%plot3(wavs,smooth(temp2(maxin2),0.1,'rloess'),ones(size(wavs))*max(max(temp))*1.18,'-','linewidth',2)
plot3(lam_obs,transmission_curves_allz(:,26),ones(size(lam_obs))*max(max(temp))*1.2,'w-','linewidth',2)
axis([3800,9000,0,1,min(min(temp)),max(max(temp))*1.5])
xlabel('\lambda [A]','fontsize',18)
ylabel('T^{IGM}_{\lambda ,z}','fontsize',18)
title('IGM transmission at z=4.5 [White line: Inoue(2014)]','fontsize',18)



%--------------------------------------------------------------------------
%----------Modeling stochastic deviations and their effects on an SED------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% getting LSST filter curves

u_crv = load('lsst_filter_curves/LSST-LSST.u.dat');
g_crv = load('lsst_filter_curves/LSST-LSST.g.dat');
r_crv = load('lsst_filter_curves/LSST-LSST.r.dat');
i_crv = load('lsst_filter_curves/LSST-LSST.i.dat');
z_crv = load('lsst_filter_curves/LSST-LSST.z.dat');
y_crv = load('lsst_filter_curves/LSST-LSST.y.dat');

filcurve(:,1) = interp1(u_crv(:,1),u_crv(:,2),wavs);
filcurve(:,2) = interp1(g_crv(:,1),g_crv(:,2),wavs);
filcurve(:,3) = interp1(r_crv(:,1),r_crv(:,2),wavs);
filcurve(:,4) = interp1(i_crv(:,1),i_crv(:,2),wavs);
filcurve(:,5) = interp1(z_crv(:,1),z_crv(:,2),wavs);
filcurve(:,6) = interp1(y_crv(:,1),y_crv(:,2),wavs);

filcurve_longer(:,1) = interp1(u_crv(:,1),u_crv(:,2),lam_obs);
filcurve_longer(:,2) = interp1(g_crv(:,1),g_crv(:,2),lam_obs);
filcurve_longer(:,3) = interp1(r_crv(:,1),r_crv(:,2),lam_obs);
filcurve_longer(:,4) = interp1(i_crv(:,1),i_crv(:,2),lam_obs);
filcurve_longer(:,5) = interp1(z_crv(:,1),z_crv(:,2),lam_obs);
filcurve_longer(:,6) = interp1(y_crv(:,1),y_crv(:,2),lam_obs);

plot(wavs,filcurve,'k');hold on
plot(lam_obs,filcurve_longer,'k:')

%--------------------------------------------------------------------------
% and a couple of galaxy SEDs to work with at z=4:

load('z4_spectra.mat')

semilogy(lam_obs,spec3_interp_lam_obs);hold on
semilogy(lam_obs,spec4_interp_lam_obs);hold on
semilogy(lam_obs,spec5_interp_lam_obs);hold on
semilogy(lam_obs,spec6_interp_lam_obs);hold on

spec2_interp_lam_obs = interp1(lam*(1+4),spec2*1e10,lam_obs);
spec2_interp_wavs = interp1(lam*(1+4),spec2*1e10,wavs);

spec3_interp_lam_obs = interp1(lam*(1+4),spec3*1e10,lam_obs);
spec3_interp_wavs = interp1(lam*(1+4),spec3*1e10,wavs);

spec4_interp_lam_obs = interp1(lam*(1+4),spec4*1e10,lam_obs);
spec4_interp_wavs = interp1(lam*(1+4),spec4*1e10,wavs);

spec5_interp_lam_obs = interp1(lam*(1+4),spec5*1e10,lam_obs);
spec5_interp_wavs = interp1(lam*(1+4),spec5*1e10,wavs);

spec6_interp_lam_obs = interp1(lam*(1+4),spec6*1e10,lam_obs);
spec6_interp_wavs = interp1(lam*(1+4),spec6*1e10,wavs);

plot(lam_obs,spec2_interp_lam_obs,'linewidth',2);hold on
plot(lam_obs,spec3_interp_lam_obs,'linewidth',2);hold on
plot(lam_obs,1e3*spec4_interp_lam_obs,'linewidth',2);hold on
plot(lam_obs,log10(1e4*spec5_interp_lam_obs),'linewidth',2);hold on
plot(lam_obs,filcurve_longer,'k:')

cvs = magma(30);

plot(lam_obs,spec4_interp_lam_obs,':');hold on
plot(wavs,spec4_interp_wavs)

load('z4_analysis.mat')
bigarr(bigarr>1) =1;
plot(wavs,bigarr(:,1))


for i = 1:20
    subplot(4,5,i)
    plot(wavs,smooth(bigarr(:,i),0.01,'loess'));hold on
    plot(wavs,smooth(bigarr(:,i+20),0.01,'loess'));hold on
    plot(wavs,smooth(bigarr(:,i+40),0.01,'loess'));hold on
    axis([3800,9000,0,1])
end

[hm1,hm2] = min(nansum(bigarr(:,1:548)))
[hm1,hm2] = max(nansum(bigarr(:,1:548)))

plot(wavs,smooth(bigarr(:,158),0.01,'loess'));hold on
plot(wavs,smooth(bigarr(:,348),0.01,'loess'));hold on

%--------------------------------------------------------------------------
figure
temp = (1e4*spec5_interp_lam_obs + spec3_interp_lam_obs);
temp2 = (1e4*spec5_interp_wavs + spec3_interp_wavs);

bgarrid = 348; % 3, 158, 348

yyaxis left
plot(lam_obs,transmission_curves_allz(:,21),'--','linewidth',1);hold on
plot(wavs,smooth(bigarr(:,bgarrid),0.01,'loess'),'color',[1,1,1]*0.9);hold on
axis([3800,9000,0,1])
xlabel('\lambda [A]','fontsize',18)
ylabel('T^{IGM}_{z,\lambda}','fontsize',18)

yyaxis right
plot(lam_obs,temp,'-','linewidth',2,'color',cvs(1,:));hold on
plot(wavs,temp2.*bigarr(:,bgarrid)','b:')
plot(lam_obs,temp.*transmission_curves_allz(:,21)','r--','linewidth',2);hold on
a = temp.*transmission_curves_allz(:,21)';
sed1 = [nansum(a.*filcurve_longer(:,2)')/nansum(filcurve_longer(:,2)),nansum(a.*filcurve_longer(:,3)')/nansum(filcurve_longer(:,3)),nansum(a.*filcurve_longer(:,4)')/nansum(filcurve_longer(:,4)),nansum(a.*filcurve_longer(:,5)')/nansum(filcurve_longer(:,5))]
b = temp2.*bigarr(:,bgarrid)';
sed2 = [nansum(b.*filcurve(:,2)')/nansum(filcurve(:,2)), nansum(b.*filcurve(:,3)')/nansum(filcurve(:,3)), nansum(b.*filcurve(:,4)')/nansum(filcurve(:,4)), nansum(b.*filcurve(:,5)')/nansum(filcurve(:,5))]
centers_lsst =[3664.37, 4807.02, 6209.82, 7542.84, 8700.52, 9633.00];
plot(centers_lsst(2:4),sed1(1:end-1),'r.','markersize',36);hold on
plot(centers_lsst(2:4),sed2(1:end-1),'b.','markersize',36)
axis([3800,9000,0,max(temp2)])
ylabel('F_\nu^{SED}','fontsize',18)

yyaxis left
plot(lam_obs,filcurve_longer/5,'k:')

%--------------------------------------------------------------------------




temp = (1e4*spec5_interp_lam_obs + spec2_interp_lam_obs);
plot(lam_obs,temp/max(temp),':','linewidth',2,'color',cvs(20,:));hold on
plot(lam_obs,temp/max(temp).*transmission_curves_allz(:,21)','-','linewidth',2,'color',cvs(20,:));hold on

plot(lam_obs,transmission_curves_allz(:,21),'k--','linewidth',2)
axis([3800,9000,0,1])


%--------------------------------------------------------------------------
% try at z=4.5

load('z4.5_spectra.mat')

semilogy(lam_obs,spec3_interp_lam_obs);hold on
semilogy(lam_obs,spec4_interp_lam_obs);hold on
semilogy(lam_obs,spec5_interp_lam_obs);hold on
semilogy(lam_obs,spec6_interp_lam_obs);hold on

spec2_interp_lam_obs = interp1(lam*(1+4.5),spec2*1e10,lam_obs);
spec2_interp_wavs = interp1(lam*(1+4.5),spec2*1e10,wavs);

spec3_interp_lam_obs = interp1(lam*(1+4.5),spec3*1e10,lam_obs);
spec3_interp_wavs = interp1(lam*(1+4.5),spec3*1e10,wavs);

spec4_interp_lam_obs = interp1(lam*(1+4.5),spec4*1e10,lam_obs);
spec4_interp_wavs = interp1(lam*(1+4.5),spec4*1e10,wavs);

spec5_interp_lam_obs = interp1(lam*(1+4.5),spec5*1e10,lam_obs);
spec5_interp_wavs = interp1(lam*(1+4.5),spec5*1e10,wavs);

spec6_interp_lam_obs = interp1(lam*(1+4.5),spec6*1e10,lam_obs);
spec6_interp_wavs = interp1(lam*(1+4.5),spec6*1e10,wavs);

plot(lam_obs,spec2_interp_lam_obs,'linewidth',2);hold on
plot(lam_obs,spec3_interp_lam_obs,'linewidth',2);hold on
plot(lam_obs,1e3*spec4_interp_lam_obs,'linewidth',2);hold on
plot(lam_obs,log10(1e4*spec5_interp_lam_obs),'linewidth',2);hold on
plot(lam_obs,filcurve_longer,'k:')

cvs = magma(30);

plot(lam_obs,spec4_interp_lam_obs,':');hold on
plot(wavs,spec4_interp_wavs)

load('z4.5_analysis.mat')
bigarr(bigarr>1) =1;
plot(wavs,bigarr(:,1))


for i = 1:20
    subplot(4,5,i)
    plot(wavs,smooth(bigarr(:,i),0.01,'loess'));hold on
    plot(wavs,smooth(bigarr(:,i+20),0.01,'loess'));hold on
    plot(wavs,smooth(bigarr(:,i+40),0.01,'loess'));hold on
    axis([3800,9000,0,1])
end

[hm1,hm2] = min(nansum(bigarr(:,1:ct_gal)))
[hm1,hm2] = max(nansum(bigarr(:,1:ct_gal)))

plot(wavs,smooth(bigarr(:,159),0.01,'loess'));hold on
plot(wavs,smooth(bigarr(:,181),0.01,'loess'));hold on

%--------------------------------------------------------------------------
figure
temp = (1e4*spec5_interp_lam_obs + spec3_interp_lam_obs);
temp2 = (1e4*spec5_interp_wavs + spec3_interp_wavs);

bgarrid = 181; % 159, 181

yyaxis left
plot(lam_obs,transmission_curves_allz(:,21),'--','linewidth',1);hold on
plot(wavs,smooth(bigarr(:,bgarrid),0.01,'loess'),'color',[1,1,1]*0.9);hold on
axis([3800,9000,0,1])
xlabel('\lambda [A]','fontsize',18)
ylabel('T^{IGM}_{z,\lambda}','fontsize',18)

yyaxis right
plot(lam_obs,temp,'-','linewidth',2,'color',cvs(1,:));hold on
plot(wavs,temp2.*bigarr(:,bgarrid)','b:')
plot(lam_obs,temp.*transmission_curves_allz(:,21)','r--','linewidth',2);hold on
a = temp.*transmission_curves_allz(:,21)';
sed1 = [nansum(a.*filcurve_longer(:,2)')/nansum(filcurve_longer(:,2)),nansum(a.*filcurve_longer(:,3)')/nansum(filcurve_longer(:,3)),nansum(a.*filcurve_longer(:,4)')/nansum(filcurve_longer(:,4)),nansum(a.*filcurve_longer(:,5)')/nansum(filcurve_longer(:,5))]
b = temp2.*bigarr(:,bgarrid)';
sed2 = [nansum(b.*filcurve(:,2)')/nansum(filcurve(:,2)), nansum(b.*filcurve(:,3)')/nansum(filcurve(:,3)), nansum(b.*filcurve(:,4)')/nansum(filcurve(:,4)), nansum(b.*filcurve(:,5)')/nansum(filcurve(:,5))]
centers_lsst =[3664.37, 4807.02, 6209.82, 7542.84, 8700.52, 9633.00];
plot(centers_lsst(2:4),sed1(1:end-1),'r.','markersize',36);hold on
plot(centers_lsst(2:4),sed2(1:end-1),'b.','markersize',36)
axis([3800,9000,0,max(temp2)])
ylabel('F_\nu^{SED}','fontsize',18)

yyaxis left
plot(lam_obs,filcurve_longer/5,'k:')

%--------------------------------------------------------------------------




temp = (1e4*spec5_interp_lam_obs + spec2_interp_lam_obs);
plot(lam_obs,temp/max(temp),':','linewidth',2,'color',cvs(20,:));hold on
plot(lam_obs,temp/max(temp).*transmission_curves_allz(:,21)','-','linewidth',2,'color',cvs(20,:));hold on

plot(lam_obs,transmission_curves_allz(:,21),'k--','linewidth',2)
axis([3800,9000,0,1])
