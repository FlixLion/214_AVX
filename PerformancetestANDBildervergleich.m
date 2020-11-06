clear all
close all

disp('#### START ####');
compStartProgram = tic;

workspaceName = datestr(now,'dd-mm-yy_HH-MM');

%Parameter
vergleichsTest = true;
performanceTest = false;
repeats = 3;
NnImg = single(4);  %Wie viele unterschiedliche Bildananazahlen
resImgMix = [10:10:200,220:20:700];

%Codes that will be compared
leonCodes = {@L0sb,@L1sb,@L2sb,@L3sb};
paulCodes = {@pL0, @pL1, @pL2, @pL3};
addsigCodes = {@addsig2vol_3_ASM, @addsig2vol_3_C};
addsigMtCodes = {@addsig2vol_3_ASM_MT, @addsig2vol_3_C_MT};
maxMt = 12;

ourCodes = [leonCodes, paulCodes];
noMtCodes= [leonCodes, paulCodes, addsigCodes];
allCodes = [leonCodes, paulCodes, addsigCodes, addsigMtCodes];
nNoMtCodes = size(noMtCodes,2);
nMtCodes = size(addsigMtCodes,2);
noMtNames = 2*size(ourCodes,2) + 2;

addsigCodes{1}(1);
addsigCodes{2}(1);
addsigMtCodes{1}(1);
addsigMtCodes{2}(1);

codeNames = string.empty;
for i=1:size(ourCodes,2)
    codeNames(i)=strcat(func2str(ourCodes{i}), "-Double");
end
currentCodeNameSize = size(codeNames,2);
for i=1:size(ourCodes,2)
    codeNames(i + currentCodeNameSize)=strcat(func2str(ourCodes{i}), "-Single");
end
currentCodeNameSize = size(codeNames,2);
for i = size(ourCodes,2)+1:size(noMtCodes,2)
    codeNames(currentCodeNameSize + i - size(ourCodes,2)) = func2str(allCodes{i});
end
currentCodeNameSize = size(codeNames,2);
for i = 0:nMtCodes - 1
    for c = 1:maxMt
        codeNames(currentCodeNameSize + i*maxMt + c) = strcat(func2str(addsigMtCodes{i+1}), " ", int2str(c),"C");
    end
end
for i = 1:size(codeNames,2)
    disp(codeNames(i));
end

codeAufzahlung = string.empty;
for i = 1:repeats
    codeAufzahlung(i) = num2str(i);
end


%Inputvariabeln
speed = single(1480);               %Geschwindigkeit des Schalls in m/s
addsigSpeed = repmat(speed, 1, NnImg);
resRl = single(1);                %BildgrÃ¶ÃŸe in Metern




if performanceTest
    times = single(zeros(repeats,size(resImgMix,2),size(codeNames,2))); %Speicher für die erzeilten Zeiten
    for i = 1:size(resImgMix,2)
        %Inputvariablen
        %Bild
        resImg = single(resImgMix(i));               %BildgrÃ¶ÃŸe in Pixeln
        pixDensity = resRl / resImg;    %Pixeldichte in Meter/Pixel
        resint = pixDensity;
        image = zeros(resImg,resImg,resImg);%Bildmatrix
        try
            imagesforaddsig32 = uint32([resImg,resImg,resImg]);
            imagesforaddsig = double(zeros(resImg,resImg,resImg));
        catch e
        end
        img_start = [0.0,0.0,0.0];
        %Scan
        posObj = [(rand()*resRl),(rand()*resRl),(rand()*resRl)];    %Random Objekt Position
        maxAscan = 2 * sqrt(sqrt(resRl^2+resRl^2)^2+resRl^2);       %m
        minAscan = 2 * sqrt(sqrt(pixDensity^2+pixDensity^2)^2+pixDensity^2);    %m
        maxTime = maxAscan / speed;         %s
        minTime = minAscan / speed;         %s
        intAscan = minTime;                 %s
        samplesAscan = 1 + maxTime / minTime;
        presetAscan = 0:intAscan:maxTime;
        speedAscan = 1 / (speed*intAscan);
        invPixDensity = 1/pixDensity;   %FÃ¼r Indexing bei Bilderstellung
        posSens = zeros(NnImg,3);
        posRecs = zeros(NnImg,3);
        Ascans = zeros(round(2 * samplesAscan),NnImg);
        for x = 1:NnImg
            posRec = [(rand()*resRl),(rand()*resRl),(rand()*resRl)];       %Random Position Reciever festelgen
            posRecs(x,:) = posRec;
            %aRecs(:,:,x) = single(posRec);
            posSen = [(rand()*resRl),(rand()*resRl),(rand()*resRl)];       %Random Position Sender festelgen
            posSens(x,:) = posSen;
            %aSens(:,:,x) = single(posSen);
            %         fprintf('%g\n',posSen)
            disGes = sqrt(sum((posSen - posObj).^2)) + sqrt(sum((posObj - posRec).^2));  %Distanzberechnung
            indexPeak = (disGes / speed) / intAscan;        %Berechnung Peak Index
            Ascan = randn(1,round(2 * samplesAscan));    %Ascan mit Rauschen kreieren
            %Ascan = zeros(1,round(2 * samplesAscan));
            Ascan(round(indexPeak)) = Ascan(round(indexPeak)) + (100/NnImg);
            Ascans(:,x) = Ascan;
        end
        addsigAscans = double(Ascans);
        
        fprintf('--- Res: %g(%g)/%g(%g) ---\n',resImgMix(i), i, resImgMix(end), size(resImgMix,2));
        for c = 1:size(codeNames,2)
            try
                fprintf(codeNames(c) + "\n");
                for n = 1:repeats
                    if c<=size(leonCodes,2)
                        tic;
                        allCodes{c}(resImg, NnImg, Ascans, speed, intAscan, posSens, posRecs, image, img_start, pixDensity);
                        times(n, i, c) = toc;
                    elseif c<=size(ourCodes,2)
                        tic;
                        allCodes{c}(resImg, NnImg, posRecs, posSens, image, Ascans, speedAscan, pixDensity, invPixDensity);
                        times(n, i, c) = toc;
                    elseif c<=size(ourCodes,2) + size(paulCodes,2)
                        tic;
                        allCodes{c-size(ourCodes,2)}(single(resImg), single(NnImg), single(Ascans), single(speed), single(intAscan), single(posSens), single(posRecs), single(image), single(img_start), single(pixDensity));
                        times(n, i, c) = toc;
                    elseif c<=2*size(ourCodes,2)
                        tic;
                        allCodes{c-size(ourCodes,2)}(single(resImg), single(NnImg), single(posRecs), single(posSens), single(image), single(Ascans), single(speedAscan), single(pixDensity), single(invPixDensity));
                        times(n, i, c) = toc;
                    elseif c<=noMtNames
                        tic;
                        allCodes{c-size(ourCodes,2)}(addsigAscans,single(img_start'),single(posRecs'),single(posSens'),single(speed),single(pixDensity),intAscan,imagesforaddsig32',imagesforaddsig);
                        times(n, i, c) = toc;
                    elseif c <= noMtNames + maxMt
                        if n==1
                            allCodes{nNoMtCodes + 1}(c-noMtNames);
                        end
                        tic;
                        allCodes{nNoMtCodes + 1}(addsigAscans,single(img_start'),posRecs',posSens',speed,pixDensity,intAscan,imagesforaddsig32',imagesforaddsig);
                        times(n, i, c) = toc;
                    elseif c <= noMtNames + 2*maxMt
                        if n==1
                            allCodes{nNoMtCodes + 2}(c-noMtNames-maxMt);
                        end
                        tic;
                        allCodes{nNoMtCodes + 2}(addsigAscans,single(img_start'),posRecs',posSens',speed,pixDensity,intAscan,imagesforaddsig32',imagesforaddsig);
                        times(n, i, c) = toc;
                    end
                    fprintf([' in ',num2str(times(n, i, c)),' Sekunden\n']);
                end
            catch e
                times(:, i, c) = NaN;
                fprintf('CRASHED\n');
                fprintf(1,'The identifier was:\n%s',e.identifier);
                fprintf(1,'There was an error! The message was:\n%s\n',e.message);
            end
        end
        save("Doku_" + workspaceName + "Unfertig");
    end
    compEndProgram = toc(compStartProgram);
    
    nCodes = size(codeNames,2);
    meanTimes = reshape(mean(times,1), [size(times,2),size(times,3)])';
    medianTimes = reshape(median(times,1), [size(times,2),size(times,3)])';
    meanDoubles = sum(meanTimes(1:8,:));
    meanSingles = sum(meanTimes(9:16,:));
    
    for i = 1:size(resImgMix,2)
        voxels(i) = NnImg*(resImgMix(i)^3);
    end
    
    CPUSpeed = 3.6e9; %GHz
    meanPerformanceClocksPerVoxel = (meanTimes./(1/CPUSpeed))./repmat(voxels,size(meanTimes,1),1);            %1./repmat(voxels,size(meanTimes,1),1)./(meanTimes./(1/CPUSpeed));
    meanPerformanceClocksPerVoxel2 = (1/CPUSpeed)./(meanTimes./repmat(voxels,size(meanTimes,1),1));
    meanPeformance = repmat(voxels,size(meanTimes,1),1)./(meanTimes.*1e9);
    
    matlabIndex = 1:nCodes-2*maxMt-2;
    AsmMt = nCodes-2*maxMt+1:nCodes-maxMt;
    CMt = nCodes-maxMt+1:nCodes;
    AsmCMt = nCodes-2*maxMt+1:nCodes;
    
    startData = 10;
    endData = 520;
    
    figure('Name','Clocks/Voxel - Mean','NumberTitle','off');
    plot(resImgMix,meanPerformanceClocksPerVoxel2(AsmMt,:));
    title('Clocks/Voxel - Mean');
    xlabel('Auflösung in Pixel');
    ylabel('Clocks/Voxel');
    legend(codeNames(AsmMt),'Location','northwest')
    
    figure('Name','Zeit - Mean','NumberTitle','off');
    imagesc(resImgMix,[1:nCodes],meanTimes);
    %imagesc(meanTimes);
    title('Zeit - Mean');
    xlabel('Auflösung in Pixel');
    ylabel('Codes');
    legend(codeNames,'Location','northwest')
    
    figure('Name','Peformance - Mean - ns','NumberTitle','off');
    imagesc(resImgMix(find(resImgMix==startData):find(resImgMix==endData)),[1:nCodes],meanPeformance(:,find(resImgMix==startData):find(resImgMix==endData)));
    title('Peformance - Mean - ns');
    xlabel('Auflösung in Pixel');
    ylabel('Codes');
    %legend(codeNames,'Location','northwest')
    
    figure('Name','Peformance - Mean - ns - ASM','NumberTitle','off');
    imagesc(resImgMix(find(resImgMix==startData):find(resImgMix==endData)),[AsmMt],meanPeformance(AsmMt,find(resImgMix==startData):find(resImgMix==endData)));
    title('Peformance - Mean - ns - ASM');
    xlabel('Auflösung in Pixel');
    ylabel('Codes');
    legend(codeNames,'Location','northwest')
    
    figure('Name','Peformance - Mean - ns - C','NumberTitle','off');
    imagesc(resImgMix(find(resImgMix==startData):find(resImgMix==endData)),[CMt],meanPeformance(CMt,find(resImgMix==startData):find(resImgMix==endData)));
    title('Peformance - Mean - ns - C');
    xlabel('Auflösung in Pixel');
    ylabel('Codes');
    legend(codeNames,'Location','northwest')
    
    figure('Name','Peformance - Mean - ns - ASM & C','NumberTitle','off');
    imagesc(resImgMix(find(resImgMix==startData):find(resImgMix==endData)),[AsmCMt],meanPeformance(AsmCMt,find(resImgMix==startData):find(resImgMix==endData)));
    title('Peformance - Mean - ns - ASM & C');
    xlabel('Auflösung in Pixel');
    ylabel('Codes');
    legend(codeNames,'Location','northwest')
    
    figure('Name','Peformance - Mean','NumberTitle','off');
    title('Peformance - Mean');
    xlabel('Auflösung in Pixel');
    ylabel('Geschwindigkeit in Voxel/Sekunde');
    hold on
    for y = 1:nCodes
        plot(resImgMix,voxels(1,:)./meanTimes(y,:))
    end
    legend(codeNames,'Location','northwest')
    hold off
    
    figure('Name','Peformance - Mean','NumberTitle','off');
    title('Peformance - Mean');
    xlabel('Auflösung in Pixel');
    ylabel('Geschwindigkeit in Voxel/Sekunde');
    plot(resImgMix,voxels./[meanSingles;meanDoubles])
    legend(["Single", "Doubles"],'Location','northwest')
    
    figure('Name','Peformance - Mean - MatlabCodes','NumberTitle','off');
    title('Peformance - Mean - MatlabCodes');
    xlabel('Auflösung in Pixel');
    ylabel('Geschwindigkeit in Voxel/Sekunde');
    hold on
    for y = matlabIndex
        plot(resImgMix,voxels(1,:)./meanTimes(y,:))
    end
    legend(codeNames(matlabIndex),'Location','northwest')
    hold off
    
    for y = 1:8
        figure('Name',['Peformance - Mean - Si vs Db - ', func2str(allCodes{y})],'NumberTitle','off');
        title(['Peformance - Mean - Si vs Db - ', func2str(allCodes{y})]);
        xlabel('Auflösung in Pixel');
        ylabel('Geschwindigkeit in Voxel/Sekunde');
        plot(resImgMix,voxels(1,:)./meanTimes([y,y+8],:))
        legend(codeNames([y,y+8]),'Location','northwest')
    end
    
    figure('Name','Peformance - Mean - ASM','NumberTitle','off');
    title('Peformance - Mean - ASM');
    xlabel('Auflösung in Pixel');
    ylabel('Geschwindigkeit in Voxel/Sekunde');
    hold on
    for y = AsmMt
        plot(resImgMix,voxels(1,:)./meanTimes(y,:))
    end
    legend(codeNames(AsmMt),'Location','northwest')
    hold off
    
    figure('Name','Peformance - Mean - C','NumberTitle','off');
    title('Peformance - Mean - C');
    xlabel('Auflösung in Pixel');
    ylabel('Geschwindigkeit in Voxel/Sekunde');
    hold on
    for y = CMt
        plot(resImgMix,voxels(1,:)./meanTimes(y,:))
    end
    legend(codeNames(CMt),'Location','northwest')
    hold off
    
    figure('Name','3D Peformance begrenzt','NumberTitle','off')
    surf(meanPeformance(:,find(resImgMix==startData):find(resImgMix==endData)));
    shading flat
    figure('Name','3D Zeiten','NumberTitle','off')
    surf(meanTimes);
    shading flat
end


if vergleichsTest
    %Inputvariablen
    %Bild
    resImg = 100;               %BildgrÃ¶ÃŸe in Pixeln
    pixDensity = resRl / resImg;    %Pixeldichte in Meter/Pixel
    resint = pixDensity;
    image = zeros(resImg,resImg,resImg);%Bildmatrix
    bilder = zeros(resImg,resImg,resImg,size(codeNames,2)); %Speicher für die erzeilten Zeiten
    try
        imagesforaddsig32 = uint32([resImg,resImg,resImg]);
        imagesforaddsig = double(zeros(resImg,resImg,resImg));
    catch e
    end
    img_start = [0.0,0.0,0.0];
    %Scan
    posObj = [(rand()*resRl),(rand()*resRl),(rand()*resRl)];    %Random Objekt Position
    maxAscan = 2 * sqrt(sqrt(resRl^2+resRl^2)^2+resRl^2);       %m
    minAscan = 2 * sqrt(sqrt(pixDensity^2+pixDensity^2)^2+pixDensity^2);    %m
    maxTime = maxAscan / speed;         %s
    minTime = minAscan / speed;         %s
    intAscan = minTime;                 %s
    samplesAscan = 1 + maxTime / minTime;
    presetAscan = 0:intAscan:maxTime;
    speedAscan = 1 / (speed*intAscan);
    invPixDensity = 1/pixDensity;   %FÃ¼r Indexing bei Bilderstellung
    posSens = zeros(NnImg,3);
    posRecs = zeros(NnImg,3);
    Ascans = zeros(round(2 * samplesAscan),NnImg);
    for x = 1:NnImg+1
        posRec = [(rand()*resRl),(rand()*resRl),(rand()*resRl)];       %Random Position Reciever festelgen
        posRecs(x,:) = posRec;
        %aRecs(:,:,x) = single(posRec);
        posSen = [(rand()*resRl),(rand()*resRl),(rand()*resRl)];       %Random Position Sender festelgen
        posSens(x,:) = posSen;
        %aSens(:,:,x) = single(posSen);
        %         fprintf('%g\n',posSen)
        disGes = sqrt(sum((posSen - posObj).^2)) + sqrt(sum((posObj - posRec).^2));  %Distanzberechnung
        indexPeak = (disGes / speed) / intAscan;        %Berechnung Peak Index
        Ascan = randn(1,round(2 * samplesAscan));    %Ascan mit Rauschen kreieren
        %Ascan = zeros(1,round(2 * samplesAscan));
        Ascan(round(indexPeak)) = Ascan(round(indexPeak)) + (100/NnImg);
        Ascans(:,x) = Ascan;
    end
    addsigAscans = double(Ascans);
    for c = 1:size(codeNames,2)
        try
            fprintf(codeNames(c) + "\n");
                if c<=size(leonCodes,2)
                    bilder(:,:,:,c)= allCodes{c}(resImg, NnImg, Ascans(:,1:NnImg), speed, intAscan, posSens(1:NnImg,:), posRecs(1:NnImg,:), image, img_start, pixDensity);
                elseif c<=size(ourCodes,2)
                    bilder(:,:,:,c)= allCodes{c}(resImg, NnImg, posRecs(1:NnImg,:), posSens(1:NnImg,:), image, Ascans(:,1:NnImg), speedAscan, pixDensity, invPixDensity);
                elseif c<=size(ourCodes,2) + size(paulCodes,2)
                    bilder(:,:,:,c)= allCodes{c-size(ourCodes,2)}(single(resImg), single(NnImg), single(Ascans(:,1:NnImg)), single(speed), single(intAscan), single(posSens(1:NnImg,:)), single(posRecs(1:NnImg,:)), single(image), single(img_start), single(pixDensity));
                elseif c<=2*size(ourCodes,2)
                    bilder(:,:,:,c)= allCodes{c-size(ourCodes,2)}(single(resImg), single(NnImg), single(posRecs(1:NnImg,:)), single(posSens(1:NnImg,:)), single(image), single(Ascans(:,1:NnImg)), single(speedAscan), single(pixDensity), single(invPixDensity));
                elseif c<=noMtNames
                    bilder(:,:,:,c)= allCodes{c-size(ourCodes,2)}(addsigAscans,single(img_start'),single(posRecs'),single(posSens'),single(speed),single(pixDensity),intAscan,imagesforaddsig32',imagesforaddsig);
                elseif c <= noMtNames + maxMt
                    allCodes{nNoMtCodes + 1}(c-noMtNames);
                    bilder(:,:,:,c)= allCodes{nNoMtCodes + 1}(addsigAscans,single(img_start'),single(posRecs'),single(posSens'),single(speed),single(pixDensity),intAscan,imagesforaddsig32',imagesforaddsig);
                elseif c <= noMtNames + 2*maxMt
                    allCodes{nNoMtCodes + 2}(c-noMtNames-maxMt);
                    bilder(:,:,:,c)= allCodes{nNoMtCodes + 2}(addsigAscans,single(img_start'),single(posRecs'),single(posSens'),single(speed),single(pixDensity),intAscan,imagesforaddsig32',imagesforaddsig);
                end
        catch e
            fprintf('CRASHED\n');
            fprintf(1,'The identifier was:\n%s',e.identifier);
            fprintf(1,'There was an error! The message was:\n%s\n',e.message);
        end
    end
    for c = 1:size(codeNames,2)
        figure('Name',codeNames(c));
        [M,I] = max(bilder(:,:,:,c));
        [M2,I2] = max(M);
        [M3,I3] = max(M2);
        imagesc([0,100],[0,100],bilder(:,:,I3,c))
        %                 bildsnr = max(data(:))/std(data(:));
        %                 erwsnr= sqrt(count)*2/mean(std(ascan));
        title(['Maximaler Wert auf ', num2str(I3*resint),' | ',num2str(I2(I3)*resint),' | ',num2str(I(I2(I3))*resint),' snr: ',num2str(M3/mean(mean(mean(bilder(:,:,:,c)))) )]);
        xlabel('x in Meter');
        ylabel('y in Meter');
    end
    c=1;
    vergleichsbild = abs(bilder(:,:,:,c)-bilder(:,:,:,c+8));
    figure('Name',[codeNames(c)+' - '+codeNames(c+8)]);
    [M,I] = max(vergleichsbild);
    [M2,I2] = max(M);
    [M3,I3] = max(M2);
    imagesc([0,100],[0,100],vergleichsbild(:,:,I3))
    %                 bildsnr = max(data(:))/std(data(:));
    %                 erwsnr= sqrt(count)*2/mean(std(ascan));
    title([codeNames(c)+' - '+codeNames(c+8)+ '  |  Maximaler Wert auf '+ num2str(I3*resint)+' | '+num2str(I2(I3)*resint)+' | '+num2str(I(I2(I3))*resint)]); %' snr: ',num2str(M3/mean(mean(mean(vergleichsbild))) )]);
    xlabel('x in Pixel');
    ylabel('y in Pixel');
end
