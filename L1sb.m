function data = L1sb(res, count, Ascans, speed, int, posSens, posRecs, data, img_start, resint)
%Enthält später die Bildinformationen

pixelpostmp = [img_start(1),img_start(2),img_start(3)];
pixelpos = repmat(pixelpostmp,res.^2,1);

multi=single(0);
for num=1:(res.^2)
    multi = multi + 1;
    pixelpos(num,3) = pixelpos(num,3)+multi.*resint;
    if multi == res
        multi = 0;
    end
end
multi = single(0);
for num=1:res.^2
    if mod(num,res) == 1
        multi = multi + 1;
    end
    pixelpos(num,2) = pixelpos(num,2)+multi.*resint;
end

for i=1:count
    ascan = Ascans(:,i);
    senderpos = posSens(i,:);
    recpos = posRecs(i,:);
    
    senderposb = repmat(senderpos,res.^2,1);
    recposb = repmat(recpos,res.^2,1);
    
    for x=1:res
        pixelpos(:,1) = img_start(1) + x.*resint;
        dist_sender_pixel = sqrt(sum((senderposb-pixelpos).^2,2));
        dist_receiver_pixel = sqrt(sum((recposb-pixelpos).^2,2));
        dges = (dist_sender_pixel + dist_receiver_pixel);
        ascanpos = round((dges/speed)/int);
        ascanpos = reshape(ascanpos,[res,res])';
        ascantmp(1,:,:) = ascan(ascanpos);
        data(x,:,:) = data(x,:,:) + ascantmp(1,:,:);
    end
    
end
% endprogram = toc(start);
% xb=[0 raumabmessungen];
% yb=[0 raumabmessungen];
% imagesc(xb,yb,data(:,:,round((objektpos(3) - img_start(3))/resint)));
% bildsnr = max(data(:))/std(data(:));
% erwsnr= sqrt(count)*2/mean(std(ascan));
% title(['erwarteter snr:', num2str(erwsnr) ,', eigentlicher snr: ' ,num2str(bildsnr)]);
% xlabel('x in m');
% ylabel('y in m');
end

%Beschriftung und Anzeigen der Messdaten
