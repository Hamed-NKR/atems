% Pair Correlation Method (PCM)
% Ramin Dastanpour & Steven N. Rogak
% Developed at the University of British Columbia
% Last updated in May 2018
% Image processing package for the analysis of TEM images. Automatic
% aggregate detection and automatic primary particle sizing


%% Clearing data and closing open windows
clear
clc;  % Clear command window
close all; % Close all figure windows except those created by imtool
imtool close all;   % Close all figure windows created by imtool


%% Choose appropriate value for based on your Excel version
xls_sheet = 2; % Uncomment this if your default MS Excel version is 2013
% xls_sheet = 4; % Uncomment this if your default MS Excel version is <2013


%% initializing values
fontSize        = 10;
minparticlesize = 4.9; % to filter out noises
% Coefficient for automatic Hough transformation
coeff_matrix    = [0.2 0.8 0.4 1.1 0.4;0.2 0.3 0.7 1.1 1.8;...
    0.3 0.8 0.5 2.2 3.5;0.1 0.8 0.4 1.1 0.5];
moreaggs    = 0;
savecounter = 0;


%% Housekeeping
% global mainfolder img img_dir FileName
extracted_text = cell(1,1);


%% Report title
report_title = {'Image_ID','Particle Area (nm^2)','Particle Perimeter (nm)',...
    'Particle Length (nm)','Particle Width (nm)','Particle aspect ratio',...
    'Radius of Gyration (nm)','Particle eq. da (nm)',...
    'dp (nm) [simple PCF]','dp (nm) [generalized PCF]','Max resolution (nm)','Engine Type'};


%% Loading images; getting image file name and directory
img.num = 0; % 0: no image loaded; 1: at least one image loaded
mainfolder = cd; % getting the directory of the code
% loop continues until at least image is selected or the program is stopped

while img.num == 0
    img_dir = '..\Images'; % get the directory of the image
    
    message = sprintf('Please choose image(s) to be analyzed');
    
    [img.files,img_dir] = uigetfile({'*.tif;*.jpg',...
        'TEM image (*.tif;*.jpg)'},'Select Images',img_dir,'MultiSelect',...
        'on');% User browses for images. Modify for other image formats
    img.num = size(img.num,2);
    if iscell(img.files) == 1 % Handling when only one image is selected
        img.files = img.files';
    elseif isempty(img.files) == 1 
        error('No image was selected');
    end
    
    if img.num == 0 % No image is selected
        
        pixsize_choise=questdlg('No image was selected! Do you want to try again?', ...
            'Error','Yes','No. Quit debugging','Yes');
        if strcmp(pixsize_choise,'No. Quit debugging')
            uiwait(msgbox('No image was selected and user decided to stop the program'))
            error('No image was selected and user decided to stop the program');
        end
    end
end
[img.num,~] = size(img.files); % Total number of images loaded


%% Main image processing loop
for img_counter = 1:img.num % run loop as many times as images selected
    
    %% Step1: Image preparation
    %% Step1-1: Loading images one-by-one
    % cd(img_dir); % change active directory to image directory
    if img.num == 1
        FileName = char(img.files); 
    else
        FileName = char(img.files(img_counter,1));
    end
    img.RawImage = imread([img_dir,FileName]); % read in image
    
    %{
    %% Step1-2: Detecting Magnification and/or pixel size
    % Determining image magnification. Image magnification can be detected
    % automatically if images are taken at UBC. In this case, image footer
    % is scanned for the defalut text and numbers used by our Quartz PCI
    % software by which TEM images are captured. Users can develop similar
    % scripts for their own images.
    % If images are not taken at UBC, and they all are taken at the same
    % magnification, this part can be modified to prevent
    % inserting/detecting magnification/pixel size for each image; and use
    % and constant value for all images.
    pixsize_choise = questdlg('Determine pixel size','Image source',...
        'Use scale bar',...
        'Insert pixel size manually',...
        'UBC -> Automatic detection');
    
    if strcmp(pixsize_choise,'Use scale bar')
        % manually choosing the magnification
        uiwait(msgbox('Please crop the image close enough to the magnification bar'))
        img.mag_crop = imcrop(img.RawImage); % crop image
        close (gcf);
        imshow(img.mag_crop); % Show Cropped image
        set(gcf,'Position',get(0,'Screensize')); % Maximize figure.
        hold on
        uiwait(msgbox('Click on a point at the start (left) of the scale bar, then on a point at the end (right) of the scale bar'));
        % user chooses two point on the edge of magnification bar
        clear bar.x bar.y
        [bar.x,bar.y] = ginput(2);
        % calculate number of pixels of magnification bar
        bar.l = abs(bar.x(2)-bar.x(1));
        line([bar.x(1),bar.x(2)],[bar.y(1),bar.y(1)],'linewidth', 3);
        dlg_title = 'Length of the magnification bar';
        promt1 = {'Please insert the length of the magnification bar in nm:'};
        num_lines = 1;
        default_l = {'100'}; %default value for user input
        % user input execution
        bar.size = str2double(cell2mat(inputdlg(promt1,dlg_title,num_lines,default_l)));
        clear num_lines dlg_title promt1 default_l
        pixsize = bar.size/bar.l;
        hold off
        close all
    elseif strcmp(pixsize_choise,'Insert pixel size manually')
        close (gcf);
        img.mag_crop = imcrop(img.RawImage); %crop image
        close (gcf);
        imshow(img.mag_crop); % Show Cropped image
        set(gcf,'Position',get(0,'Screensize')); % Maximize figure.
        dlg_title = 'Pixel size';
        promt1 = {'Please insert the pixel size in nm/pixel:'};
        num_lines = 1;
        default_l = {'0.535592'}; %default value for user input
        % user input execution
        pixsize = str2double(cell2mat(inputdlg(promt1,dlg_title,num_lines,default_l)));
    end
    %}
    
    %% Step 1-3: Crop footer and get scale from footer
    
    [img,pixsize] = tools.get_scale_img(img);
    
    % Build the image processing coefficients for the image based on its
    % magnification
    if pixsize <= 0.181
        coeffs = coeff_matrix(1,:);
    elseif pixsize <= 0.361
        coeffs = coeff_matrix(2,:);
    else 
        coeffs = coeff_matrix(3,:);
    end
    % displaying the image
    imshow(img.RawImage);
    title('processing Image', 'FontSize', fontSize);
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
    
    
    %% Step 1-4: Saving Cropped image
    close (gcf);
    
    %checking whether the Output folder is available 
    dirName = sprintf('Data/PCMOutput');
    
    if exist(dirName,'dir') ~= 7 % 7 if exist parameter is a directory
        mkdir(dirName) % make output folder
    end
    
%     % save Cropped image as a .tif file. Change to desired format if needed.
%     imwrite...
%     (img.Cropped,[FileName '_CroppedImage_' num2str(img_counter) '.tif'])
%     cd(mainfolder) % return to the code directry
    

    %% Step 2: automatic/semi-automatic aggregate detection
    img.Binary = PCM.Agg_detection(img,pixsize,moreaggs,minparticlesize,coeffs);
    
    img.Edge        = edge(img.Binary,'sobel'); % Sobel edge detection
    SE              = strel('disk',1);
    img.DilatedEdge = imdilate(img.Edge,SE); % morphological dilation
    
    clear img.Edge SE2
    img.Imposed = imimposemin(img.Cropped, img.DilatedEdge);
    
    
    %% Step 3: Automatic primary particle sizing
    %% Step 3-1: Find and size all particles on the final binary image
    CC = bwconncomp(abs(img.Binary-1));
    NofAggs = CC.NumObjects; % count number of particles
    
    %% Lopp for each aggregate on the image under investigation
    for nAgg = 1:NofAggs
        %% Step 3-2: Prepare an image of the isolated aggregate
        Agg.Image = zeros(size(img.Binary));
        Agg.Image(CC.PixelIdxList{1,nAgg}) = 1;
        
        % Edge Image via Sobel
        % Use Sobel's Method as a built-in edge detection function for
        % aggregates's outline. Other methods (e.g. Roberts, Canny)can also
        % be considered
        Image_edge = edge(Agg.Image,'sobel');
        
        % Dilated Edge Image
        % Use dilation to strengthen the aggregate's outline
        SE = strel('disk',1);
        Dilated_Image_edge = imdilate(Image_edge,SE);
        clear Image_edge SE
        FinalImposedImage = imimposemin(img.Cropped, Dilated_Image_edge);
        figure
        hold on
        imshow(FinalImposedImage);
%         %save Binary image as a .tif file. Can be changed to other formats
%         imwrite(Agg.Image,[FileName '_BinaryImage_' num2str(img_counter) '.tif'])
        
        %% Step 3-3: Development of the pair correlation function (PCF)
        
        %% 3-3-1: Find the skeleton of the aggregate
        Skel = bwmorph(Agg.Image,'thin',Inf);
        [skeletonY, skeletonX] = find(Skel);
        
        %% 3-3-2: Calculate the distances between skeleton pixels and other pixels
        clear row col
        [row, col] = find (Agg.Image);
        
        % Total number of pixels in aggregate
        Agg.Pixels = nnz(Agg.Image);
        
        % to consolidate the pixels of consideration to much smaller arrays
        
        p       = 0;
        X       = 0;
        Y       = 0;
        density = 20; % larger densities make the program less computationally expensive
        
        for i = 1:density:Agg.Pixels
            p    = p + 1;
            X(p) = col(i);
            Y(p) = row(i);
        end
        
        % to calculate all the distances with reference to one pixel at a time,
        % using vector algebra, and then adding the results
        for j = 1:1:length(skeletonX)
            Distance_int = ((X-skeletonX(j)).^2+(Y-skeletonY(j)).^2).^.5;
            Distance_mat(((j-1)*length(X)+1):(j*length(X))) = Distance_int;
        end
        
        %% 3-3-3: Construct the pair correlation
        % Sort radii into bins and calculate PCF
        Distance_max = double(uint16(max(Distance_mat)));
        Distance_mat = nonzeros(Distance_mat).*pixsize;
        nbins        = Distance_max * pixsize;
        dr           = 1;
        Radius       = 1:dr:nbins;
        
        % Pair correlation function (PCF)
        PCF = hist(Distance_mat, Radius);
        
        % Smoothing the pair correlation function (PCF)
        d                                  = 5 + 2*Distance_max;
        BW1                                = zeros(d,d);
        BW1(Distance_max+3,Distance_max+3) = 1;
        BW2                                = bwdist(BW1,'euclidean');
        BW4                                = BW2./Distance_max;
        BW5                                = im2bw(BW4,1);
        
        clear row col
        
        [row,col]            = find(~BW5);
        Distance_Denaminator = ((row-Distance_max+3).^2+(col-Distance_max+3).^2).^.5;
        Distance_Denaminator = nonzeros(Distance_Denaminator).*pixsize;
        Denamonator          = hist(Distance_Denaminator,Radius);
        Denamonator          = Denamonator.*length(skeletonX)./density;
        PCF                  = PCF./Denamonator;
        PCF_smoothed         = smooth(PCF);
        
        for i=1:size(PCF_smoothed)-1
            if PCF_smoothed(i) == PCF_smoothed(i+1)
                PCF_smoothed(i+1) = PCF_smoothed(i)-1e-4;
            end
        end
        
        %% 3-4: Computing aggregate dimentions/parameters
        % Projected area of the aggregate
        Agg.Area = Agg.Pixels*pixsize^2;
        
        % Projected area equivalent diameter of the aggregate
        Agg.da = ((Agg.Pixels/pi)^.5)*2*pixsize;
        
        % Aggregate Perimeter
        clear row col
        [row, col]           = find(Agg.Image);
        perimeter_pixelcount = 0;
        
        for k = 1:1:Agg.Pixels
            if      (Agg.Image(row(k)-1,col(k))==0) || ...
                    (Agg.Image(row(k),col(k)+1)==0) || ...
                    (Agg.Image(row(k),col(k)-1)==0) || ...
                    (Agg.Image(row(k)+1,col(k))==0)
                
                perimeter_pixelcount = perimeter_pixelcount + 1;
            end
        end
        
        Agg.Perimeter = perimeter_pixelcount*pixsize;
        
        % Aggregate gyration radius
        [xpos,ypos] = find(Agg.Image);
        n_pix       = length(xpos);
        Centroid.x  = sum(xpos)/n_pix;
        Centroid.y  = sum(ypos)/n_pix;
        Ar          = zeros(n_pix,1);
        
        for k = 1:n_pix;
            Ar(k,1) = (((xpos(k,1)-Centroid.x)*pixsize)^2+...
                      ((ypos(k,1)-Centroid.y)*pixsize)^2)*pixsize^2;
        end
        
        Agg.Rg = (sum(Ar)/Agg.Area)^0.5;
        clear Ar ypos xpos Centroid.x Centroid.y
        
        % Aggregate Length, Width, and aspect ratio      
        Agg.L           = max((max(row)-min(row)),(max(col)-min(col)))*pixsize;
        Agg.W           = min((max(row)-min(row)),(max(col)-min(col)))*pixsize;
        Agg.Aspectratio = Agg.L/Agg.W;
        
        %% 3-5: Primary particle sizing
        %% 3-5-1: Simple PCM
        PCF_simple   = .913;
        Agg.RpSimple = interp1(PCF_smoothed, Radius, PCF_simple);
        
        %% 3-5-2: Generalized PCM
        URg      = 1.1*Agg.Rg; % 10% higher than Rg
        LRg      = 0.9*Agg.Rg; % 10% lower than Rg
        PCFRg    = interp1(Radius, PCF_smoothed, Agg.Rg); % P at Rg
        PCFURg   = interp1(Radius, PCF_smoothed, URg); % P at URg
        PCFLRg   = interp1(Radius, PCF_smoothed, LRg); % P at LRg
        PRgslope = (PCFURg+PCFLRg-PCFRg)/(URg-LRg); % dp/dr(Rg)
        
        PCF_generalized   = (.913/.84)*(0.7+0.003*PRgslope^(-0.24)+0.2*Agg.Aspectratio^-1.13);
        Agg.RpGeneralized = interp1(PCF_smoothed, Radius, PCF_generalized);
        
        %% Plot pair correlation function in line graph format
        filename = 'Pair Correlation Function Plot.jpeg';
        path     = sprintf('PCF_%d_%d_new',nAgg,img_counter);
        str      = sprintf('Pair Correlation Line Plot %f ',PCF_simple);
        
        figure, loglog(Radius, smooth(PCF), '-r'), title (str), xlabel ('Radius'), ylabel('PCF(r)')
        
        hold on
        
        loglog(Agg.RpSimple,PCF_simple,'*')
   %    saveas(gcf, path) %uncomment to save
        close all
        
        %% Clear variables
        clear Agg.Image Image_edge FinalImposedImage Skel skeletonY ...
            skeletonX Agg.Pixels Distance_mat Radius PCF PCF_smoothed ...
            Denamonator row col
        
        
        %% Step 4: Saving results
        close all
        extracted_data(1) = round(Agg.Area,2);
        extracted_data(2) = round(Agg.Perimeter,2);
        extracted_data(3) = round(Agg.L,2);
        extracted_data(4) = round(Agg.W,2);
        extracted_data(5) = round(Agg.Aspectratio,2);
        extracted_data(6) = round(Agg.Rg,2);
        extracted_data(7) = round(Agg.da,2);
        extracted_data(8) = round(2*Agg.RpSimple,2); % dp from simple PCM
        extracted_data(9) = round(2*Agg.RpGeneralized,2); % dp from generalized PCM
        extracted_data(10) = round(pixsize,3);
        extracted_text(1)={FileName};

        %% Step 4-1: Autobackup
        if exist('Data\PCMOutput\PCM_data.mat','file')==2
            save('Data\PCMOutput\PCM_data.mat','extracted_data','extracted_text','-append');
        else
            save('Data\PCMOutput\PCM_data.mat','extracted_data','extracted_text','report_title');
        end
        if exist('Data\PCMOutput\PCM_Output.xls','file')==2
            [~,sheets,~] = xlsfinfo('Data\PCMOutput\PCM_Output.xls'); %finding info about the excel file
            sheetname=char(sheets(1,xls_sheet)); % choosing the second sheet
            [datanum ~]=xlsread('Data\PCMOutput\PCM_Output.xls',sheetname); %loading the data
            starting_row=size(datanum,1)+2;
            xlswrite('Data\PCMOutput\PCM_Output.xls',extracted_text,'TEM_Results',['A' num2str(starting_row)]);
            xlswrite('Data\PCMOutput\PCM_Output.xls',extracted_data,'TEM_Results',['B' num2str(starting_row)]);
        else
            savecounter=1;
            xlswrite('Data\PCMOutput\PCM_Output.xls',report_title,'TEM_Results','A1');
            xlswrite('Data\PCMOutput\PCM_Output.xls',extracted_text,'TEM_Results','A2');
            xlswrite('Data\PCMOutput\PCM_Output.xls',extracted_data,'TEM_Results','B2');
        end
        
    end
    
    
end
