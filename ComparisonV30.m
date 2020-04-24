%V20 übergibt mehr Argumente

clear all
close all
%Programm

disp('#### START ####');
compStartProgram = tic;

workspaceName =  datestr(now,'mmmm-dd-yyyy_HH.MM.SS.FFF');

repeats = 1;
NnImg = 3;  %Wie viele unterschiedliche Bildananazahlen

NResImg = 50;  %Wie viele unterschiedliche Bildauflösungen
NResImg2 = 35;  %Wie viele unterschiedliche Bildauflösungen
NResImg3 = 20;
stepResImg = 10; %resImg Stepsize
stepResImg2 = 20; %resImg Stepsize
stepResImg3 = 50;
resImgMix = [stepResImg:stepResImg:NResImg * stepResImg, NResImg * stepResImg + stepResImg2:stepResImg2:NResImg2 * stepResImg2, NResImg2 * stepResImg2 + stepResImg3:stepResImg3:700];
resImgMix = (100);


addsigCode = {@addsig2vol_3};
leonCodes = {@L0sb,@L1sb,@L2sb,@L3sb};
paulCodes = {@pL0, @pL1, @pL2, @pL3};
allCodes = [leonCodes, paulCodes, addsigCode];
nCodes = size(leonCodes,2) + size(paulCodes,2) + size(addsigCode,2);
codeNames = string.empty;
codeAufzahlung = string.empty;
for i = 1:size(allCodes,2)
    codeNames(i) = func2str(allCodes{i});
end
for i = 1:repeats
    codeAufzahlung(i) = num2str(i);
end

times = zeros(repeats,size(resImgMix,2),nCodes); %Speicher für die erzeilten Zeiten
voxels = times(:,:,1);

for i = 1:size(resImgMix,2)
    for n = 1:repeats
        voxels(n,i) = NnImg*(resImgMix(i)^3);
    end
end


for i = 1:size(resImgMix,2)
    images = (zeros(resImgMix(i),resImgMix(i),resImgMix(i),nCodes));
    speed = (1480);               %Geschwindigkeit des Schalls in m/s
    
    %Bild
    resImg = (resImgMix(i));               %BildgrÃ¶ÃŸe in Pixeln
    resRl =(1.5);                %BildgrÃ¶ÃŸe in Metern
    pixDensity = (resRl / resImg);    %Pixeldichte in Meter/Pixel
    resint = pixDensity;
    image = zeros(resImg,resImg,resImg);    %Bildmatrix
    
    img_start = ([0.0,0.0,0.0]);
    
    %Scan
    rng('default');
    posObj = [(rand()*resRl),(rand()*resRl),(rand()*resRl)];    %Random Objekt Position
    maxAscan = 2 * sqrt(sqrt(resRl^2+resRl^2)^2+resRl^2);       %m
    minAscan = 2 * sqrt(sqrt(pixDensity^2+pixDensity^2)^2+pixDensity^2);    %m
    maxTime = maxAscan / speed;         %s
    minTime = minAscan / speed;         %s
    intAscan = minTime;                 %s
    samplesAscan = maxTime / minTime;
    presetAscan = 0:intAscan:maxTime;
    speedAscan = 1 / (speed*intAscan);
    invPixDensity = 1/pixDensity;   %FÃ¼r Indexing bei Bilderstellung
    
    posSens = (zeros(NnImg,3));
    posRecs = (zeros(NnImg,3));
    Ascans = zeros(round(samplesAscan),NnImg);
    for x = 1:NnImg
        posRec = ([(rand()*resRl),(rand()*resRl),(rand()*resRl)]);       %Random Position Reciever festelgen
        posRecs(x,:) = (posRec);
        posSen = ([(rand()*resRl),(rand()*resRl),(rand()*resRl)]);       %Random Position Sender festelgen
        posSens(x,:) = (posSen);
        %fprintf('%g\n',posSen)
        disGes = sqrt(sum((posSen - posObj).^2)) + sqrt(sum((posObj - posRec).^2));  %Distanzberechnung
        indexPeak = (disGes / speed) / intAscan;        %Berechnung Peak Index
        Ascan = randn(1,round(samplesAscan));    %Ascan mit Rauschen kreieren
        %Ascan = zeros(1,round(samplesAscan));
        Ascan(round(indexPeak)) = Ascan(round(indexPeak)) + (100/NnImg);
        Ascans(:,x) = Ascan;
    end
    
    x=resImg;

    fprintf('--- Res: %g(%g)/%g(%g) ---\n',resImgMix(i), i, resImgMix(end), size(resImgMix,2));
    for c = 1
        for n = 1:repeats
            for Nimg = 1:NnImg
%                 try
                    fprintf(func2str(addsigCode{c}));
                    tic;
                    images(:,:,:,c+size(leonCodes,2) + size(paulCodes,2)) = addsigCode{c}(Ascans(:,Nimg),single(reshape(img_start,3,1)),single(reshape(posRecs(Nimg,:),3,1)),single(reshape(posSens(Nimg,:),3,1)),single(speed),single(resImg),single(intAscan),uint32([x,x,x]),images(:,:,:,c+nCodes));
%                     images(:,:,:,c+size(leonCodes,2) + size(paulCodes,2)) = zeros(resImgMix(i),resImgMix(i),resImgMix(i)); %addsigCode{c}(Ascans(:,Nimg),single(reshape(img_start,3,1)),single(reshape(posRecs(Nimg,:),3,1)),single(reshape(posSens(Nimg,:),3,1)),single(speed),single(resImg),single(intAscan),uint32([x,x,x]),images(:,:,:,c+nCodes));
                    
                    times(n,i,c+size(leonCodes,2) + size(paulCodes,2)) = toc;
                    fprintf([' in ',num2str(times(n,i,c+size(leonCodes,2) + size(paulCodes,2))),' Sekunden\n']);
%                 catch e
%                     times(n,i,c+size(leonCodes,2) + size(paulCodes,2)) = NaN;
%                     fprintf(' CRASHED\n');
%                     fprintf(1,'The identifier was:\n%s',e.identifier);
%                     fprintf(1,'\n%s',e.message);
%                 end
            end
        end
    end
    
    for c = 1:size(leonCodes,2)
        for n = 1:repeats
            try
                fprintf(func2str(leonCodes{c}));
                tic;
                images(:,:,:,c) = leonCodes{c}(resImg, NnImg, Ascans, speed, intAscan, posSens, posRecs, image, img_start, pixDensity);
                times(n,i,c) = toc;
                fprintf([' in ',num2str(times(n,i,c)),' Sekunden\n']);
            catch e
                times(n,i,c) = NaN;
                fprintf(' CRASHED\n');
                fprintf(1,'The identifier was:\n%s',e.identifier);
                fprintf(1,'\n%s',e.message);
            end
        end
    end
    
    
    for c = 1:size(paulCodes,2)
        for n = 1:repeats
            try
                fprintf(func2str(paulCodes{c}));
                tic;
                
                images(:,:,:,c+size(leonCodes,2)) = paulCodes{c}(resImg, NnImg, posRecs, posSens, image, Ascans, speedAscan, pixDensity, invPixDensity);
                times(n,i,c+size(leonCodes,2)) = toc;
                fprintf([' in ',num2str(times(n,i,c+size(leonCodes,2))),' Sekunden\n']);
            catch e
                times(n,i,c+size(leonCodes,2)) = NaN;
                fprintf(' CRASHED\n');
                fprintf(1,'The identifier was:\n%s',e.identifier);
                fprintf(1,'\n%s',e.message);
            end
        end
    end
    save(['CompWspaceUnfertig_',workspaceName,'.mat'])
end
for c = 1:size(allCodes,2)
    figure('Name',func2str(allCodes{c}));
    [M,I] = max(images(:,:,:,c));
    [M2,I2] = max(M);
    [M3,I3] = max(M2);
    imagesc([0,resImgMix(end)],[0,resImgMix(end)],images(:,:,I3,c))
    %                 bildsnr = max(data(:))/std(data(:));
    %                 erwsnr= sqrt(count)*2/mean(std(ascan));
    title(['Maximaler Wert auf ', num2str(I3*resint),' | ',num2str(I2(I3)*resint),' | ',num2str(I(I2(I3))*resint),' snr: ',num2str(M3/mean(mean(mean(images(:,:,:,c)))) )]);
    xlabel('x in Meter');
    ylabel('y in Meter');
end

maxfehler = [nCodes,nCodes,2];
for i = 1:nCodes
    for j = 1:nCodes
        vergleich = abs(images(:,:,:,i) - images(:,:,:,j));
        %hist(vergleich(:),1000);
        [M,I] = max(vergleich);
        [M2,I2] = max(M);
        [M3,I3] = max(M2);
        maxfehler(i,j,1:2) = [M3,I3];
        if i < j
            figure('Name',['Bilder Vergleich: ',func2str(allCodes{i}),' vs ', func2str(allCodes{j})]);
            imagesc([0,resImgMix(end)],[0,resImgMix(end)],vergleich(:,:,I3))
            title(['Maximaler Wert auf ', num2str(I3*resint),' | ',num2str(I2(I3)*resint),' | ',num2str(I(I2(I3))*resint),' snr: ',num2str(M3/mean(vergleich(:)))]);
            xlabel('x in Meter');
            ylabel('y in Meter');
        end
    end
end
figure('Name','Bilder Vergleich Übersicht');
imagesc(maxfehler(:,:,1))
meanTimes = mean(times,1);
medianTimes = median(times,1);

compEndProgram = toc(compStartProgram);

%disp(times)
%close all
figure('Name','Geschwindigkeitsvergleich - 1');
title('Geschwindigkeitsvergleich - 1');
xlabel('Auflösung in Pixel');
ylabel('Zeit in Sekunde');
hold on
for y = 1:size(paulCodes,2)+size(leonCodes,2)     %Zahl auf Anzahl der Codes änder !!!!!!
    plot(resImgMix,meanTimes(1,:,y))
end
legend(codeNames,'Location','northwest')
hold off

figure('Name','Peformance - Mean - log');
set(gca, 'YScale', 'log');
title('Peformance - Mean');
xlabel('Auflösung in Pixel');
ylabel('Geschwindigkeit in Voxel/Sekunde');
hold on
for y = 1:size(paulCodes,2)+size(leonCodes,2)     %Zahl auf Anzahl der Codes änder !!!!!!
    semilogy(resImgMix,voxels(1,:)./meanTimes(1,:,y))
end
legend(codeNames,'Location','northwest')
hold off

figure('Name','Peformance - Median - log');
set(gca, 'YScale', 'log');
title('Peformance - Median');
xlabel('Auflösung in Pixel');
ylabel('Geschwindigkeit in Voxel/Sekunde');
hold on
for y = 1:size(paulCodes,2)+size(leonCodes,2)     %Zahl auf Anzahl der Codes änder !!!!!!
    semilogy(resImgMix,voxels(1,:)./medianTimes(1,:,y))
end
legend(codeNames,'Location','northwest')
hold off

figure('Name','Peformance - Mean');
title('Peformance - Mean');
xlabel('Auflösung in Pixel');
ylabel('Geschwindigkeit in Voxel/Sekunde');
hold on
for y = 1:size(paulCodes,2)+size(leonCodes,2)     %Zahl auf Anzahl der Codes änder !!!!!!
    plot(resImgMix,voxels(1,:)./meanTimes(1,:,y))
end
legend(codeNames,'Location','northwest')
hold off

figure('Name','Peformance - Median');
title('Peformance - Median');
xlabel('Auflösung in Pixel');
ylabel('Geschwindigkeit in Voxel/Sekunde');
hold on
for y = 1:size(paulCodes,2)+size(leonCodes,2)     %Zahl auf Anzahl der Codes änder !!!!!!
    plot(resImgMix,voxels(1,:)./medianTimes(1,:,y))
end
legend(codeNames,'Location','northwest')
hold off

figure('Name','Peformance - Durchgänge Vergleich');
title('Peformance - Durchgänge Vergleich');
xlabel('Auflösung in Pixel');
ylabel('Geschwindigkeit in Voxel/Sekunde');
hold on
for y = 1:repeats
    plot(resImgMix,voxels(1,:)./mean(times(y,:,:),3))
end
legend(codeAufzahlung,'Location','northwest')
hold off

fprintf('\n#### Alles in %g Sekunden (%g Stunden) ####',compEndProgram,compEndProgram/60/60)

save(['CompWspace_',workspaceName,'.mat'])