function [debug,image] = pL2(resImg, nImg, posSens, posRecs, image, Ascans, speedAscan, pixDensity, invPixDensity, debugOn)
debug = NaN;
if debugOn == true
    debug = (zeros(resImg,resImg,resImg));
end
line = (zeros(resImg,3));     %Eine Linie
line(:,3) = pixDensity:pixDensity:resImg*pixDensity;    %Letzte Koordinate durchloopen

for i = 1:nImg
    Ascan = Ascans(:,i);
    posSen = posSens(i,:);
    posRec = posRecs(i,:);
    for x = pixDensity:pixDensity:resImg*pixDensity
        line(:,1) = repmat(x,resImg,1);
        xI = round((x * invPixDensity));
        
        for y = pixDensity:pixDensity:resImg*pixDensity
            line(:,2) = repmat(y,resImg,1);
            yI = round((y*invPixDensity));
            lineDisGes = sqrt(sum((repmat(posSen,resImg,1)-line).^2,2)) + sqrt(sum((line-repmat(posRec,resImg,1)).^2,2));   %Berechnung Sen-Pix(Line) + Pix(Line)-Rec
            if debugOn == true
                debug(xI,yI,:) = lineDisGes;
            end
            image(xI,yI,:) = image(xI,yI,:) + reshape(Ascan(round(lineDisGes*speedAscan))',[1,1,resImg]);                       %Bild überschreiben
        end
    end
end
end