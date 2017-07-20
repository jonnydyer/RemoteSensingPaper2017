% Make some figures for HDR section of Jonny's paper "Imaging Trends and the Future 
% of Satellite Remote Sensing"
%
% 2017-07-18, mark robertson

clear

symbols = {':', '-.', '-'};

do_gifs = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNR curve when combining 16 images

% Range of digital numbers
dn = [1:2:4096];

% Standard deviation of noise follows sigma_scale * sqrt(dn)
sigma_scale = 0.7;
snr_single = dn./(sigma_scale*sqrt(dn));

figure(1), clf
num_combine = 16;
plot(dn, snr_single), xlabel('Digital Number (DN)'),ylabel('SNR'),hold on
snr_many = sqrt(num_combine)*snr_single;
plot(dn,snr_many)
title('SNR, Photon Shot Noise Dominated'),grid off
legend({'SNR Single Image',['SNR ' num2str(num_combine) ' Images']},'Location','Best')
axis tight

set(gcf,'Units','inches')
p = get(gcf,'Position');
set(gcf,'Position',[p(1),p(2),8,6])
set(gcf, 'PaperPositionMode','auto')
print('-dpng', '-r600',['hdr-design-snr-' num2str(num_combine) '-images.png']);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNR curve when considering different exposures

figure(2),clf

% Scale for variance, instead of for standard deviation
sigma2_scale = sigma_scale^2;

gains = [1,2,4];

dn_norm = dn/max(dn);
leg={};
sigma2num=zeros(size(dn));
sigma2den=zeros(size(dn));
for index=1:length(gains)
  gain = gains(index)
  last = floor(length(dn)/gain);
  plot(dn_norm(1:last),sqrt(gain)*snr_single(1:last),symbols{index}),hold on
  leg{end+1} = [num2str(gain) 'x Exposure'];
  beta = zeros(size(dn)); 
  beta(1:last) = gain;
  sigma2num = sigma2num + beta.*sigma2_scale.*dn;
  sigma2den = sigma2den + beta;
end
sigma = sqrt(sigma2num)./sigma2den;
plot(dn_norm,dn./sigma,'k--')
leg{end+1} = 'Net SNR After Combination';
ylim([0,150])
legend(leg,'Location','NorthEast')
xlabel('Normalized Brightness to 1x')
ylabel('SNR')
title('Compare SNR and Dynamic Range for Different Exposures')
grid off

set(gcf,'Units','inches')
p = get(gcf,'Position');
set(gcf,'Position',[p(1),p(2),8,6])
set(gcf, 'PaperPositionMode','auto')
print('-dpng', '-r600',['hdr-design-snr-hdr-examples.png']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNR curve when building up a set of equal exposures

N = 16

if do_gifs
  values = 1:N;
  for value=values
    figure(3),clf
    plot(dn_norm,sqrt(value)*snr_single),hold on
    plot(dn_norm,snr_single,'--')
    plot(dn_norm,sqrt(values(end))*snr_single,'--')
    legend({sprintf('Combining %2d Images',value),'Single Frame SNR',sprintf('Combining %2d Images',values(end))},'Location','SouthEast');
    ylim([0,300])
    xlabel('Normalized Brightness to 1x')
    ylabel('SNR')
    title('Building Up Multi-Frame SNR');
    grid off
    set(gcf,'Units','inches')
    p = get(gcf,'Position');
    set(gcf,'Position',[p(1),p(2),8,6])
    set(gcf, 'PaperPositionMode','auto')
    print('-dpng', '-r144',sprintf('hdr-temp-%02d.png',value));
  end
  system(sprintf('convert hdr-temp-??.png -delay 25 -loop 0 hdr-design-build-simple-snr.gif'));
  system(sprintf('rm hdr-temp-??.png'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Animated GIF for SNR curve when considering different exposures



%gains = [1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.026 1.131 ...
%         1.245 1.376 1.517 1.673 1.847 2.038 2.244 2.479 2.726 3.002 3.310 3.650 ...
%         4.023 4.475 4.910 5.438 5.987 6.659 7.339 7.998];
gains = [1.000 1.000 1.000 1.000 1.000 1.047 1.268 1.530 ...
         1.857 2.244 2.705 3.279 3.976 4.840 5.884 7.048];

gains = gains(1:N)       
       
sigma2num=zeros(size(dn));
sigma2den=zeros(size(dn));
% First get final SNR curve
for gain=gains
  last = floor(length(dn)/gain);
  beta = zeros(size(dn)); 
  beta(1:last) = gain;
  sigma2num = sigma2num + beta.*sigma2_scale.*dn;
  sigma2den = sigma2den + beta;
end
sigma_shaped = sqrt(sigma2num)./sigma2den;

if do_gifs
  leg={};
  sigma2num=zeros(size(dn));
  sigma2den=zeros(size(dn));
  counter=0;
  for gain=gains
    figure(4),clf
    last = floor(length(dn)/gain);
    beta = zeros(size(dn)); 
    beta(1:last) = gain;
    sigma2num = sigma2num + beta.*sigma2_scale.*dn;
    sigma2den = sigma2den + beta;

    sigma = sqrt(sigma2num)./sigma2den;
    plot(dn_norm,dn./sigma),hold on
    plot(dn_norm,snr_single,'--');
    plot(dn_norm,dn./sigma_shaped,'--')
    legend({sprintf('Combining %2d Exopsures',counter+1),'Single Frame SNR',sprintf('Combinging %2d Exposures',length(gains))},'Location','SouthEast');
    ylim([0,300])
    xlabel('Normalized Brightness to 1x')
    ylabel('SNR')
    title('Equalized SNR')
    grid off

    set(gcf,'Units','inches')
    p = get(gcf,'Position');
    set(gcf,'Position',[p(1),p(2),8,6])
    set(gcf, 'PaperPositionMode','auto')
    print('-dpng', '-r144',sprintf('hdr-temp2-%02d.png',counter));

    counter = counter+1;
  end
  system(sprintf('convert hdr-temp2-??.png -delay 25 -loop 0 hdr-design-build-hdr-snr.gif'));
  system(sprintf('rm hdr-temp2-??.png'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare standard vs. HDR curves

figure(5),clf,hold on
N = length(gains);
plot(dn_norm, snr_single, ':')
plot(dn_norm,sqrt(N)*snr_single,'--')
plot(dn_norm,dn./sigma_shaped, 'k')
legend({'Single Image','Fixed Integration Time','HDR Exposures'},'Location','SouthEast')
xlabel('Normalized Brightness to 1x')
ylabel('SNR')
title(['SNR For Standard vs. HDR Exposures, Combining ' num2str(N) ' Images'])
grid off

set(gcf,'Units','inches')
p = get(gcf,'Position');
set(gcf,'Position',[p(1),p(2),8,6])
set(gcf, 'PaperPositionMode','auto')
print('-dpng', '-r600',['hdr-design-compare-snr.png']);


