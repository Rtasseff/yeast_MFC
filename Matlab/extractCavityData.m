% Gilles Charvin - 01/2014
% Extract parameters from microfludics RLS data in order to compare with
% Ryan's simulation

% Principles : first load a project position

%TO DO :

%0) Mean cell size as a function of time DONE
%1) average density as a function of time DONE
%2) cell mean squared displacement as a function of time DONE
%3) assume bidding direction is defined by a unit vector from the cell center to the bud neck, this direction changes and in the rough 2D approximation this change can be quantified as angular displacement in the the x,y plane.  It would be great to have the mean squared angular displacement as a function of time.
%4) difusión profile (I will get this from the paper first, but am including it here for documentation purposes)


% sample project :

% load project

%path='/Users/charvin/Documents/Work/Data/movies/120405_CP03_full_spectum/120405_CP03_full_spectrum_hi-pH-project.mat';
path='/Users/charvin/Documents/Work/Data/movies/20120803/120803_mito_pH5-project.mat-project.mat';

global timeLapse
out=at_load(path);

%% load semgnetation variable for position 8;

global segmentation
pathseg=['/Users/charvin/Documents/Work/Data/movies/20120803/120803_mito_pH5-project.mat-pos8/segmentation-batch.mat'];
out=load(pathseg);

%% load first image to determine the contour of the cavity

img=phy_loadTimeLapseImage(position,1,1,'nonretreat'); % load image 1 channel phase contrast

figure, imshow(img,[]);

BW = roipoly;

close
totalcavitypixels=sum(sum(BW));

%% plot readouts size and density

global segmentation

position=segmentation.position;

totarea=zeros(1,size(segmentation.cells1,1));

meanarea=zeros(1,size(segmentation.cells1,1));

for i=1:size(segmentation.cells1,1)
    cc=0;
    for j=1:size(segmentation.cells1,2)
        if segmentation.cells1(i,j).n~=0
            totarea(i)=totarea(i)+segmentation.cells1(i,j).area;
            meanarea(i)=meanarea(i)+segmentation.cells1(i,j).area;
            cc=cc+1;
        end
        
    end
    meanarea(i)=meanarea(i)/cc;
end

totarea=totarea/totalcavitypixels;

tim=1:1:size(segmentation.cells1,1);
tim=tim*10;

figure, plot(tim,totarea); xlabel('time (min)'); ylabel('Cell density'); title('Evolution of cell density in chip as a function of time');

saveas(gca,['Celldensity-pos' num2str(position) '.fig']);
myExportFig(['Celldensity-pos' num2str(position) '.pdf']);
save(['Celldensity-pos' num2str(position) '.mat'],'totarea');

figure, plot(tim,meanarea); xlabel('time (min)'); ylabel('Cell area (pixels)'); title('Evolution of mean cell area as a function of time');

saveas(gca,['Cellsize-pos' num2str(position) '.fig']);
myExportFig(['Cellsize-pos' num2str(position) '.pdf']);
save(['Cellsize-pos' num2str(position) '.mat'],'meanarea');

%% compute Cell MSD

MSD=zeros(numel(segmentation.tcells1),1000); errMSD=zeros(1,1000); MSDavg=zeros(1,1000);
MSDang=zeros(numel(segmentation.tcells1),1000); errMSDang=zeros(1,1000); MSDavgang=zeros(1,1000);

cc2=1;
for j=1:numel(segmentation.tcells1)
    
    %if j~=698
    %    continue
    %end
    
    im=[segmentation.tcells1(j).Obj];
    [im ix]=sort(im);
    
    if length(im)>10
        
        ox=[segmentation.tcells1(j).Obj(ix).ox];
        oy=[segmentation.tcells1(j).Obj(ix).oy];
        
        
        for i=1:length(im)
           x= segmentation.tcells1(j).Obj(ix(i)).x;
           y= segmentation.tcells1(j).Obj(ix(i)).y;
           e = fitEllipse(x,y);
           ang(i)= -e.phi;
        end
        
        for i=1:length(im)/2
            m=floor(length(im)/i);
            cc=1;
            for k=1:m
                %cc2
                MSD(j,i)=MSD(j,i)+(ox(i*(k-1)+i)-ox(i*(k-1)+1)).^2+ (oy(i*(k-1)+i)-oy(i*(k-1)+1)).^2 ;
                
                MSDang(j,i)=MSDang(j,i)+ (ang(i*(k-1)+i)-ang(i*(k-1)+1)).^2 ;
                cc=cc+1;
            end
            
            %`a=MSD(cc2,1:20)
            
            % cc2
            MSD(j,i)=MSD(j,i)/(cc-1);
            MSDang(j,i)=MSDang(j,i)/(cc-1);
        end
        
        j
    end
    
end

for i=1:size(MSD,2)
    avg=MSD(:,i);
    avgang=MSDang(:,i);
    
    pix=find(avg~=0);
    
    if numel(pix)
        
    avg=avg(pix);
    MSDavg(i)=mean(avg);
    errMSD(i)=std(avg)./sqrt(length(avg));
    
    avgang=avgang(pix);
    
    MSDavgang(i)=mean(avgang);
    errMSDang(i)=std(avgang)./sqrt(length(avgang));
    
    end
end


tim=1:1:size(segmentation.cells1,1);
tim=tim*10;

lastt=find(MSDavg~=0,1,'last');
%MSDavg

figure, errorbar(tim(1:lastt),MSDavg(1:lastt),errMSD(1:lastt),'LineWidth',1); %errorbar(MSDavg,errMSD);

xlabel('time (min)'); ylabel('MSD (pixel^2)'); title('Mean square displacement of cells');

saveas(gca,['CellMSD-pos' num2str(position) '.fig']);
myExportFig(['CellMSD-pos' num2str(position) '.pdf']);
outMSD=MSDavg(1:lastt);
save(['CellMSD-pos' num2str(position) '.mat'],'outMSD');

%lastt=find(MSDavg~=0,1,'last')

figure, errorbar(tim(1:lastt),MSDavgang(1:lastt),errMSDang(1:lastt),'LineWidth',1); %errorbar(MSDavg,errMSD);
xlabel('time (min)'); ylabel('MSD (radians^2)'); title('Mean square angular displacement of cells');

saveas(gca,['CellMSDang-pos' num2str(position) '.fig']);
myExportFig(['CellMSDang-pos' num2str(position) '.pdf']);
outMSD=MSDavgang(1:lastt);
save(['CellMSD-pos' num2str(position) '.mat'],'outMSD');
