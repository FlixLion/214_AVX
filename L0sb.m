%Benchmarken: Unterschied: gleiche Input/Output Variable / anderer Name

function data = L0sb(res, count, Ascans, speed, int, posSens, posRecs, data, img_start, resint)

pixelpos = single(zeros(res,res,res,3));
pixelpos(:,:,:,1) = img_start(1);
pixelpos(:,:,:,2) = img_start(2);
pixelpos(:,:,:,3) = img_start(3);

for i = 1:res
    x = img_start(1) + i.*resint;
    pixelpos(i,:,:,1) = x;
end
for i = 1:res
    y = img_start(2) + i.*resint;
    pixelpos(:,i,:,2) = y;
end
for i = 1:res
    z = img_start(3) + i.*resint;
    pixelpos(:,:,i,3) = z;
end



for i=1:count
    
    
    ascan = Ascans(:,i);
    senderpos = posSens(i,:);
    recpos = posRecs(i,:);
    
    senderposb = single(zeros(res,res,res,3));
    senderposb(:,:,:,1) = senderpos(1);
    senderposb(:,:,:,2) = senderpos(2);
    senderposb(:,:,:,3) = senderpos(3);
    recposb = single(zeros(res,res,res,3));
    recposb(:,:,:,1) = recpos(1);
    recposb(:,:,:,2) = recpos(2);
    recposb(:,:,:,3) = recpos(3);
    
    
    dist_sender_pixel = sqrt(sum((senderposb-pixelpos).^2,4));
    dist_receiver_pixel = sqrt(sum((recposb-pixelpos).^2,4));
    dges = (dist_sender_pixel + dist_receiver_pixel);
    ascanpos = round((dges/speed)/int);
    data = data + ascan(ascanpos);
    
end
end
%
%  xb=[0 raumabmessungen];
%  yb=[0 raumabmessungen];
%  imagesc(xb,yb,data(:,:,round((objektpos(3) - img_start(3))/resint)));
%  bildsnr = max(data(:))/std(data(:));
%  erwsnr= sqrt(count)*2/mean(std(ascan));
%  title(['erwarteter snr:', num2str(erwsnr) ,', eigentlicher snr: ' ,num2str(bildsnr)]);
%  xlabel('x in m');
%  ylabel('y in m');
%Beschriftung und Anzeigen der Messdaten
