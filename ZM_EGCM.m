% function [palette]=ZM_EGCM(file,path)
% input the filename and filepath
% output the ZM amplitude and EGCM
[file, path]=uigetfile({...
'*.jpg;*.bmp;*.png;*.jpeg','All supported format(*.jpg;*.bmp;*.png;*.jpeg)';...
'*.bmp','Image format bmp (*.bmp)';...          % file selemean_seg_pixelstion with multiple
'*.jpg','Image format jpg (*.jpg)';...
'*.png','Image format png (*.png)';...
'*.jpeg','Image format jpeg (*.jpeg)';...
},'MultiSelect','on');
if ~path                                  % si rien n'est selemean_seg_pixelstionner
        return;
end
if ~iscell(file)
        temp = file;
        file = cell(1);
        file{1} = temp;
end
    NumOfFile=size(file,2);
    palette=cell(NumOfFile,1);
    steps=NumOfFile;
    step=1;
    l=1;
    order=25;
    Num_of_ZM=36;
    zm=1;
    A=zeros(Num_of_ZM,1);
    hwait=waitbar(0,'Waiting>>>>>>>>');
    while l<=NumOfFile
       Input_im = imread([path file{l}]);
       if size(Input_im,3)>2
           p=double(imresize(rgb2gray(Input_im),[256,256]));
       else
           p=double(imresize(Input_im,[256,256]));
       end
% im_re=cell(size(C,1),1);
N = size(p,1);
M = size(p,2);
recon2=zeros(N,M);
for O=0:order
    ZM_coff=(-O):2:O;
    repetition=ZM_coff(:,ZM_coff>=0);
    for Rep=1:length(repetition)
%% GET ZERNIKE MOMENT 
x = 1:M; y = 1:N;
[X,Y] = meshgrid(x,y);
Radius = sqrt((2.*(X+0.5)-N-1).^2+(2.*(Y+0.5)-N-1).^2)/(sqrt(2)*N);
Theta = atan2((N-1-2.*Y+2),(2.*X-N+1-2));
Radius = (Radius<=1).*Radius;
Rad = radialpoly(Radius,O,repetition(Rep));    % get the radial polynomial

Product = p.*Rad.*exp(-1i*repetition(Rep)*Theta);
Z = sum(Product(:));        % calculate the moments

cnt = (pi*(N-1).^2)/2;             % count the number of pixels inside the unit circle
Z = (O+1)*Z/cnt;            % normalize the amplitude of moments
A(zm,1) = abs(Z);                 % calculate the amplitude of the moment
% Phi = angle(Z)*180/pi;      % calculate the phase of the mement (in degrees)
    zm=zm+1;
    recon1=Z.*Rad.*exp(-1i*repetition(Rep)*Theta);  %% image reconstraction 
    recon2=recon2+recon1;
    end
end
zm=1;
%% EGCM
BW=edge(uint8(p),'canny');
[Gmag,Gdir] = imgradient(BW);
[Phi_map, index] = imquantize(Gdir,[-157.5,-112.5,-67.6,-22.5,22.5,67.5,112.5,157.5],[5,6,7,8,1,2,3,4,5]);  %%[-180,-135,-90,-45,0,45,90,135,-180] label with one of eight orientations to update Phi
Phi_map(BW==0)=NaN;  %% label the pixel (not the coordinates of the boundary pixels) as no gradient with value NaN.
offset=[0 1;-1 1;-1 0;-1 -1;0 -1;1 -1;1 0;1 1]; %% distance of one pixel in the direction of the eight gradient orientations.
glcm = graycomatrix(Phi_map,'Offset',offset,'Symmetric', false,'NumLevels',8,'G',[]);
EGCM=sum(glcm,3);
EGCM_nor=reshape(EGCM,[64,1]); %% each row in CM is then concatenated into a vector of length 64
EGCM_nor=(EGCM_nor-min(EGCM_nor(:)))./(max(EGCM_nor(:))-min(EGCM_nor(:))); %% and values are then normaliezed in the range of [0 1]
palette{l,1}=A;
palette{l,2}=EGCM_nor;
num_per=(step/(steps))*100;
       num_per = sprintf('%8.2f',num_per);
       str=['extracting features... ',num2str(num_per),'% Completed'];
waitbar((step/(steps)),hwait,str);                                       % updating the progress bar                                                
        pause(0.01);
        step=step+1; 
        l=l+1;
    end                                                          % incrementing progress bar step
    close(hwait) ;
% end