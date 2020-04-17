function image = pL3(resImg, nImg, posSens, posRecs, image, Ascans, speedAscan, pixDensity, useless)
for i = 1:nImg
    Ascan = Ascans(:,i);
    posSen = posSens(i,:);
    posRec = posRecs(i,:);
    posPix=(zeros(1,3));
    for x = 1:resImg
        posPix(1)=x;
        for y = 1:resImg
            posPix(2)=y;
            for z = 1:resImg
                posPix(3)=z;
                disGes = sqrt(sum((posSen-posPix.*pixDensity).^2)) + sqrt(sum((posPix.*pixDensity-posRec).^2));
                image(x,y,z) = image(x,y,z) + Ascan(round(disGes * speedAscan));
            end
        end
    end
end

end