function data = L2sb(res, count, Ascans, speed, int, posSens, posRecs, data, img_start, resint)
pixelposz = single(zeros(res,1));
for num=1:res
    pixelposz(num) = num*resint;
end
pixelpostmp = [img_start(1),img_start(2),img_start(3)];
pixelpos = repmat(pixelpostmp,res,1);
pixelpos(:,3) = pixelposz + pixelpos(:,3);

for i=1:count
    ascan = Ascans(:,i);
    senderpos = posSens(i,:);
    recpos = posRecs(i,:);
    %fprintf('%d ', i);
    
    senderposb = repmat(senderpos,res,1);
    recposb = repmat(recpos,res,1);
    
    for x=1:res
        pixelpos(:,1) = repmat((img_start(1) + x.*resint),res,1);
        for y=1:res
            pixelpos(:,2) = repmat((img_start(1) + y.*resint),res,1);
%             dist_sender_pixel = sqrt((senderposb(:,1)-pixelpos(:,1)).^2+(senderposb(:,2)-pixelpos(:,2)).^2+(senderposb(:,3)-pixelpos(:,3)).^2);
%             dist_receiver_pixel = sqrt((recposb(:,1)-pixelpos(:,1)).^2+(recposb(:,2)-pixelpos(:,2)).^2+(recposb(:,3)-pixelpos(:,3)).^2);
            dist_sender_pixel = sqrt(sum((senderposb-pixelpos).^2,2));
            dist_receiver_pixel = sqrt(sum((recposb-pixelpos).^2,2));
            dges = (dist_sender_pixel + dist_receiver_pixel);
            ascanpos = round((dges/speed)/int);
            ascantmp(1,1,:) = ascan(ascanpos);
            data(x,y,:) = data(x,y,:) + ascantmp(1,1,:);
        end
    end
    
end
end

% xb=[0 raumabmessungen];
% yb=[0 raumabmessungen];
% imagesc(xb,yb,data(:,:,round((objektpos(3) - img_start(3))/resint)));
% bildsnr = max(data(:))/std(data(:));
% erwsnr= sqrt(count)*2/mean(std(ascan));
% title(['erwarteter snr:', num2str(erwsnr) ,', eigentlicher snr: ' ,num2str(bildsnr)]);
% xlabel('x in m');
% ylabel('y in m');
% Beschriftung und Anzeigen der Messdaten
