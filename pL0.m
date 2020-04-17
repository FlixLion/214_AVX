function [debug,image] = pL0(resImg, nImg, posSens, posRecs, image, Ascans, speedAscan, pixDensity, useless, debugOn)
debug = NaN;
if debugOn == true
    debug = (zeros(resImg,resImg,resImg));
end

cube = zeros(resImg,resImg,resImg,3);    %Eine Cube mit Koordinaten
cube(:,:,:,1) = repmat([pixDensity:pixDensity:(resImg)*pixDensity]', [1, resImg, resImg]);
cube(:,:,:,2) = repmat(pixDensity:pixDensity:(resImg)*pixDensity, [resImg, 1, resImg]);
cube(:,:,:,3) = repmat(reshape([pixDensity:pixDensity:(resImg)*pixDensity],[1,1,resImg]), [resImg, resImg, 1]);

for i = 1:nImg
    Ascan = Ascans(:,i);
    posSen = posSens(i,:);
    posRec = posRecs(i,:);
    posSen3D = cat(4,posSen(1),posSen(2),posSen(3));
    posRec3D = cat(4,posRec(1),posRec(2),posRec(3));
    
    cubeDisGes = sqrt(sum(((posSen3D(ones(resImg,1),ones(resImg,1),ones(resImg,1),:))-cube).^2,4)) + sqrt(sum((cube - (posRec3D(ones(resImg,1),ones(resImg,1),ones(resImg,1),:))).^2,4));
    if debugOn == true
        debug = cubeDisGes;
    end
    image = image + Ascan(round(cubeDisGes*speedAscan));
end
end

