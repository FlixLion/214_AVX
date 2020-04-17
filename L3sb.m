function disges = L3sbDebug(res, count, Ascans, speed, int, posSens, posRecs, data, img_start, resint)
disges = single(zeros(res,res,res,count));
for i=1:count
    ascan = Ascans(:,i);
    senderpos = posSens(i,:);
    recpos = posRecs(i,:);
    pixelpos = single(zeros(3,1));
    for z=1:res
        pixelpos(3) = img_start(3) + z.*resint;
        for x=1:res
            pixelpos(1) = img_start(1) + x.*resint;
            for y=1:res
                pixelpos(2) = img_start(2) + y.*resint;
                dist_sender_pixel = sqrt((senderpos(1)-pixelpos(1)).^2+(senderpos(2)-pixelpos(2)).^2+(senderpos(3)-pixelpos(3)).^2);
                dist_receiver_pixel = sqrt((recpos(1)-pixelpos(1)).^2+(recpos(2)-pixelpos(2)).^2+(recpos(3)-pixelpos(3)).^2);
                dges = (dist_sender_pixel + dist_receiver_pixel);
                ascanpos = round((dges/speed)/int);
                data(x,y,z) = data(x,y,z) + ascan(ascanpos); %Berechnung eines der Einzelbilder und zusammenfügen mit den vorherigen
                %Einzelbildern
                disges(x,y,z,i) = dges;
            end
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
%Beschriftung und Anzeigen der Messdaten