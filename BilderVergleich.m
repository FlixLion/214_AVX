%V20 übergibt mehr Argumente
%mex -setup Compiler installieren
%disges Ausgeben und Vergleichen
%addsig2vol3 in Comparison integrieren
%-> gleiches Bild, Geschwindigkeit überprüfen
%Flags von addsig2vol3 ausprobieren (Multithreading, Assembler vs C,...)
%(Zeile 60-70)
clear all
close all
%Programm
disp('#### START ####');
compStartProgram = tic;
%rng('default');
workspaceName = datestr(now,'dd_mm_yy_HH_MM');

repeats = 1;
NnImg = single(20);  %Wie viele unterschiedliche Bildananazahlen

NResImg = 50;  %Wie viele unterschiedliche Bildauflösungen
NResImg2 = 35;  %Wie viele unterschiedliche Bildauflösungen
NResImg3 = 20;
stepResImg = 10; %resImg Stepsize
stepResImg2 = 20; %resImg Stepsize
stepResImg3 = 50;
resImgMix = [stepResImg:stepResImg:NResImg * stepResImg, NResImg * stepResImg + stepResImg2:stepResImg2:NResImg2 * stepResImg2, NResImg2 * stepResImg2 + stepResImg3:stepResImg3:700];
resImgMix = single([200]);

leonCodes = {@L0sb,@L1sb,@L2sb,@L3sb};
paulCodes = {@pL0, @pL1, @pL2, @pL3};
leonCodes = {};
paulCodes = {};
%addsigCodes = {@addsig2vol_3_C, @addsig2vol_3_ASM, @addsig2vol_3_C_MT, @addsig2vol_3_ASM_MT};
addsigCodes = {@addsig2vol_3};

addsigCodes{1}(4);
%addsigCodes{2}(1);
% addsigCodes{3}(4);
% addsigCodes{4}(4);


nCodes = size(leonCodes,2) + size(paulCodes,2) + size(addsigCodes,2);
codeNames = string.empty;
codeAufzahlung = string.empty;
for i = 1:size(leonCodes,2)
    codeNames(i) = func2str(leonCodes{i});
end
for i = size(leonCodes,2)+1:size(paulCodes,2)+size(leonCodes,2)
    codeNames(i) = func2str(paulCodes{i-size(leonCodes,2)});
end
for i = size(leonCodes,2)+1+size(paulCodes,2):size(paulCodes,2)+size(leonCodes,2)+size(addsigCodes,2)
    codeNames(i) = func2str(addsigCodes{i-size(leonCodes,2)-size(paulCodes,2)});
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
images = single(zeros(resImgMix(i),resImgMix(i),resImgMix(i),nCodes));

for i = 1:size(resImgMix,2)
    
    
    
    
    speed = single(1480);               %Geschwindigkeit des Schalls in m/s
    
    %Bild
    resImg = single(resImgMix(i));               %BildgrÃ¶ÃŸe in Pixeln
    resRl =single(1);                %BildgrÃ¶ÃŸe in Metern
    pixDensity = single(resRl / resImg);    %Pixeldichte in Meter/Pixel
    resint = pixDensity;
    image = zeros(resImg,resImg,resImg);%Bildmatrix
    imagesforaddsig32 = uint32([resImgMix(i),resImgMix(i),resImgMix(i)]);
    imagesforaddsig = double(zeros(resImgMix(i),resImgMix(i),resImgMix(i)));
    img_start = single([0.0,0.0,0.0]);

    %Scan
    posObj = [(rand()*resRl),(rand()*resRl),(rand()*resRl)];    %Random Objekt Position
    maxAscan = 2 * sqrt(sqrt(resRl^2+resRl^2)^2+resRl^2);       %m
    minAscan = 2 * sqrt(sqrt(pixDensity^2+pixDensity^2)^2+pixDensity^2);    %m
    maxTime = maxAscan / speed;         %s
    minTime = minAscan / speed;         %s
    intAscan = minTime;                 %s
    %samplesAscan = 1 + maxTime / minTime;
    samplesAscan = single(3000);
    presetAscan = 0:intAscan:maxTime;
    speedAscan = 1 / (speed*intAscan);
    invPixDensity = 1/pixDensity;   %FÃ¼r Indexing bei Bilderstellung
    
    posSens = single(zeros(NnImg,3));
    posRecs = single(zeros(NnImg,3));
    Ascans = single(zeros(round(2 * samplesAscan),NnImg));
    for x = 1:NnImg
        posRec = single([(rand()*resRl),(rand()*resRl),(rand()*resRl)]);       %Random Position Reciever festelgen
        posRecs(x,:) = single(posRec);
        posSen = single([(rand()*resRl),(rand()*resRl),(rand()*resRl)]);       %Random Position Sender festelgen
        posSens(x,:) = single(posSen);
        fprintf('%g\n',posSen)
        disGes = sqrt(sum((posSen - posObj).^2)) + sqrt(sum((posObj - posRec).^2));  %Distanzberechnung
        indexPeak = (disGes / speed) / intAscan;        %Berechnung Peak Index
        Ascan = randn(1,round(2 * samplesAscan));    %Ascan mit Rauschen kreieren
        %Ascan = zeros(1,round(2 * samplesAscan));
        Ascan(round(indexPeak)) = Ascan(round(indexPeak)) + (100/NnImg);
        Ascans(:,x) = Ascan;
    end
    addsigAscans = double(Ascans);
    
    fprintf('--- Res: %g(%g)/%g(%g) ---\n',resImgMix(i), i, resImgMix(end), size(resImgMix,2));
    for c = 1:size(leonCodes,2)
        for n = 1:repeats
            try
                fprintf(func2str(leonCodes{c}));
                tic;
                images(:,:,:,c) = leonCodes{c}(resImg, NnImg, Ascans, speed, intAscan, posSens, posRecs, image, img_start, pixDensity);
                times(n,i,c) = toc;
                fprintf([' in ',num2str(times(n,i,c)),' Sekunden\n']);
                figure('Name',func2str(leonCodes{c}));
                [M,I] = max(images(:,:,:,c));
                [M2,I2] = max(M);
                [M3,I3] = max(M2);
                imagesc([0,resImgMix(i)],[0,resImgMix(i)],images(:,:,I3,c))
                title(['Maximaler Wert auf ', num2str(I3*resint),' | ',num2str(I2(I3)*resint),' | ',num2str(I(I2(I3))*resint),' snr: ',num2str(M3/mean(mean(mean(images(:,:,:,c)))) )]);
                xlabel('x in Meter');
                ylabel('y in Meter');
            catch e
                times(n,i,c) = NaN;
                fprintf(' CRASHED\n');
                fprintf(1,'The identifier was:\n%s',e.identifier);
                fprintf(1,'There was an error! The message was:\n%s\n',e.message);
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
                figure('Name',func2str(paulCodes{c}));
                [M,I] = max(images(:,:,:,c+size(leonCodes,2)));
                [M2,I2] = max(M);
                [M3,I3] = max(M2);
                imagesc([0,resImgMix(i)],[0,resImgMix(i)],images(:,:,I3,c+size(leonCodes,2)))
                title(['Maximaler Wert auf ', num2str(I3*resint),' | ',num2str(I2(I3)*resint),' | ',num2str(I(I2(I3))*resint),' snr: ',num2str(M3/mean(mean(mean(images(:,:,:,c+size(leonCodes,2))))) )]);
                xlabel('x in Meter');
                ylabel('y in Meter');
            catch e
                times(n,i,c+size(leonCodes,2)) = NaN;
                fprintf(' CRASHED\n');
                fprintf(1,'The identifier was:\n%s',e.identifier);
                fprintf(1,'There was an error! The message was:\n%s\n',e.message);
            end
        end
    end
    for c = 1:size(addsigCodes,2)
        for n = 1:repeats
            try
                fprintf(func2str(addsigCodes{c}));
                tic;
                for b = 1:NnImg
                    images(:,:,:,c+size(leonCodes,2)+size(paulCodes,2)) = images(:,:,:,c+size(leonCodes,2)+size(paulCodes,2)) + addsigCodes{c}(addsigAscans(:,b),single(img_start'),posRecs(b,:)',posSens(b,:)',speed,pixDensity,intAscan,imagesforaddsig32',imagesforaddsig);
                end
                %image = addsig2vol_3_ASM_MT(Ascans(:,b),single(img_start'),posRecs(b,:)',posSens(b,:)',speed,pixDensity,intAscan,imagesforaddsig32',imagesforaddsig);
                times(n,i,c+size(leonCodes,2)+size(paulCodes,2)) = toc;
                fprintf([' in ',num2str(times(n,i,c+size(leonCodes,2)+size(paulCodes,2))),' Sekunden\n']);
                figure('Name',func2str(addsigCodes{c}));
                [M,I] = max(images(:,:,:,c+size(leonCodes,2)+size(paulCodes,2)));
                [M2,I2] = max(M);
                [M3,I3] = max(M2);
                imagesc([0,resImgMix(i)],[0,resImgMix(i)],images(:,:,I3,c+size(leonCodes,2)+size(paulCodes,2)))
                title(['Maximaler Wert auf ', num2str(I3*resint),' | ',num2str(I2(I3)*resint),' | ',num2str(I(I2(I3))*resint),' snr: ',num2str(M3/mean(mean(mean(images(:,:,:,c+size(leonCodes,2)+size(paulCodes,2))))) )]);
                xlabel('x in Meter');
                ylabel('y in Meter');
            catch e
                times(n,i,c) = NaN;
                fprintf(' CRASHED\n');
                fprintf(1,'The identifier was:\n%s',e.identifier);
                fprintf(1,'There was an error! The message was:\n%s\n',e.message);
            end
        end
    end
    save('ComparisonWorkspaceUnfertig_02-04')
end
figure('Name','Bilder Vergleich','NumberTitle','off');
for i = 1:nCodes
    for j = 1:nCodes
        vergleich = abs(images(:,:,:,i) - images(:,:,:,j));
        [M,I] = max(vergleich);
        [M2,I2] = max(M);
        [M3,I3] = max(M2);
        maxfehler(i,j) = M3;
        imagesc([0,100],[0,100],vergleich(:,:,I3))
        title(['Maximaler Wert auf ', num2str(I3*resint),' | ',num2str(I2(I3)*resint),' | ',num2str(I(I2(I3))*resint),' snr: ',num2str(M3/mean(vergleich(:)))]);
        xlabel('x in Meter');
        ylabel('y in Meter');
    end
end
figure('Name','Bilder Vergleich Übersicht','NumberTitle','off');
imagesc(maxfehler)
meanTimes = mean(times,1);
medianTimes = median(times,1);

compEndProgram = toc(compStartProgram);

%disp(times)
%close all
% figure('Name','Geschwindigkeitsvergleich - 1','NumberTitle','off');
% title('Geschwindigkeitsvergleich - 1');
% xlabel('Auflösung in Pixel');
% ylabel('Zeit in Sekunde');
% hold on
% for y = 1:size(paulCodes,2)+size(leonCodes,2)     %Zahl auf Anzahl der Codes änder !!!!!!
%     plot(resImgMix,meanTimes(1,:,y))
% end
% legend(codeNames,'Location','northwest')
% hold off
% 
% figure('Name','Peformance - Mean - log','NumberTitle','off');
% set(gca, 'YScale', 'log');
% title('Peformance - Mean');
% xlabel('Auflösung in Pixel');
% ylabel('Geschwindigkeit in Voxel/Sekunde');
% hold on
% for y = 1:size(paulCodes,2)+size(leonCodes,2)     %Zahl auf Anzahl der Codes änder !!!!!!
%     semilogy(resImgMix,voxels(1,:)./meanTimes(1,:,y))
% end
% legend(codeNames,'Location','northwest')
% hold off
% 
% figure('Name','Peformance - Median - log','NumberTitle','off');
% set(gca, 'YScale', 'log');
% title('Peformance - Median');
% xlabel('Auflösung in Pixel');
% ylabel('Geschwindigkeit in Voxel/Sekunde');
% hold on
% for y = 1:size(paulCodes,2)+size(leonCodes,2)     %Zahl auf Anzahl der Codes änder !!!!!!
%     semilogy(resImgMix,voxels(1,:)./medianTimes(1,:,y))
% end
% legend(codeNames,'Location','northwest')
% hold off
% 
% figure('Name','Peformance - Mean','NumberTitle','off');
% title('Peformance - Mean');
% xlabel('Auflösung in Pixel');
% ylabel('Geschwindigkeit in Voxel/Sekunde');
% hold on
% for y = 1:size(paulCodes,2)+size(leonCodes,2)     %Zahl auf Anzahl der Codes änder !!!!!!
%     plot(resImgMix,voxels(1,:)./meanTimes(1,:,y))
% end
% legend(codeNames,'Location','northwest')
% hold off
% 
% figure('Name','Peformance - Median','NumberTitle','off');
% title('Peformance - Median');
% xlabel('Auflösung in Pixel');
% ylabel('Geschwindigkeit in Voxel/Sekunde');
% hold on
% for y = 1:size(paulCodes,2)+size(leonCodes,2)     %Zahl auf Anzahl der Codes änder !!!!!!
%     plot(resImgMix,voxels(1,:)./medianTimes(1,:,y))
% end
% legend(codeNames,'Location','northwest')
% hold off
% 
% figure('Name','Peformance - Durchgänge Vergleich','NumberTitle','off');
% title('Peformance - Durchgänge Vergleich');
% xlabel('Auflösung in Pixel');
% ylabel('Geschwindigkeit in Voxel/Sekunde');
% hold on
% for y = 1:repeats
%     plot(resImgMix,voxels(1,:)./mean(times(y,:,:),3))
% end
% legend(codeAufzahlung,'Location','northwest')
% hold off
% 
% fprintf('\n#### Alles in %g Sekunden (%g Stunden) ####',compEndProgram,compEndProgram/60/60)
% 
% save(workspaceName)