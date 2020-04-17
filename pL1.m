function [debug,image] = pL1(resImg, nImg, posSens, posRecs, image, Ascans, speedAscan, pixDensity, invPixDensity, debugOn)
debug = NaN;
if debugOn == true
    debug = (zeros(resImg,resImg,resImg));
end
square = (zeros(resImg,resImg,3));    %Eine Square mit Koordinaten
for x = 1:resImg
    square(:,x,3) = pixDensity:pixDensity:(resImg)*pixDensity;
    square(:,x,2) = pixDensity * (x);
end


for i = 1:nImg
    Ascan = Ascans(:,i);
    posSen = posSens(i,:);
    posRec = posRecs(i,:);
    for x = pixDensity:pixDensity:(resImg)*pixDensity
        square(:,:,1) = repmat(x,resImg,resImg);
        xI = round((x * invPixDensity));
        squareDisGes = sqrt(sum( (repmat( cat(3,posSen(1),posSen(2),posSen(3)),resImg,resImg )-square).^2,3 )) + sqrt(sum((square - repmat(cat(3,posRec(1),posRec(2),posRec(3)),resImg,resImg)).^2,3));
        if debugOn == true
            debug(xI,:,:) = squareDisGes';
        end
        image(xI,:,:) = image(xI,:,:) + reshape(Ascan(round(squareDisGes*speedAscan,0))',[1,resImg,resImg]);
    end
end
end