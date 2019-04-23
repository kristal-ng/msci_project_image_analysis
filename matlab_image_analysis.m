% Initialize workspace
clear;
clc;


% Define size of cell you want to normalize to. This is for cells that are not the reference.
norm = 55989

% Read HDF5 files that are output from STORM analysis

hdfA = 'ANTXR.hdf5'
A.x = h5read(hdfA, '/tracks/tracks_0/x') + 1
A.y = h5read(hdfA, '/tracks/tracks_0/y') + 1
ANTXR = [A.x A.y]

hdfN = 'NUP.hdf5'
N.x = h5read(hdfN, '/tracks/tracks_0/x') + 1
N.y = h5read(hdfN, '/tracks/tracks_0/y') + 1
NUP155 = [N.x N.y]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fiducials aligned according to reference hyb0

hdfM = 'MYC.hdf5'
M.x = h5read(hdfM, '/tracks/tracks_0/x') + 1
M.y = h5read(hdfM, '/tracks/tracks_0/y') + 1
preMYC = [M.x M.y]

hdfP = 'PGRMC.hdf5'
P.x = h5read(hdfP, '/tracks/tracks_0/x') + 1
P.y = h5read(hdfP, '/tracks/tracks_0/y') + 1
prePGRMC = [P.x P.y]

%Align coordinates according to affine transform matrix for hyb1

ta = 5.4699; tb = 0.9995; tc = 0.0014;
td = 8.0649; te = 0.0011; tf = 0.9991;

T1 = [tb te 0;
 tc tf 0;
 ta td 1];

tform1 = maketform('affine', T1);

MYC = tformfwd(preMYC, tform1)
PGRMC = tformfwd(prePGRMC, tform1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hdfRV = 'RAVER.hdf5'
RV.x = h5read(hdfRV, '/tracks/tracks_0/x') + 1
RV.y = h5read(hdfRV, '/tracks/tracks_0/y') + 1
preRAVER = [RV.x RV.y]

hdfRB = 'RBFOX.hdf5'
RB.x = h5read(hdfRB, '/tracks/tracks_0/x') + 1
RB.y = h5read(hdfRB, '/tracks/tracks_0/y') + 1
preRBFOX = [RB.x RB.y]

%Align coordinates according to affine transform matrix for hyb2

tg = 8.4719; 
th = 0.9968;
ti = -0.0023;
tj = 11.3170;
tk = 0.0005;
tl = 1.0006;

T2 = [th tk 0;
      ti tl 0;
      tg tj 1];

tform2 = maketform('affine', T2);

RAVER = tformfwd(preRAVER, tform2)
RBFOX = tformfwd(preRBFOX, tform2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hdfT = 'TARDBP.hdf5'
T.x = h5read(hdfT, '/tracks/tracks_0/x') + 1
T.y = h5read(hdfT, '/tracks/tracks_0/y') + 1
preTARDBP = [T.x T.y]

hdfV = 'VDAC.hdf5'
V.x = h5read(hdfV, '/tracks/tracks_0/x') + 1
V.y = h5read(hdfV, '/tracks/tracks_0/y') + 1
preVDAC = [V.x V.y]

%Align coordinates according to affine transform matrix for hyb3

tm = 8.6751;
tn = 0.9991;
to = 0.0001;
tp = 12.6326;
tq = -0.0022;
tr = 0.9992;

T3 = [tn tq 0;
      to tr 0;
      tm tp 1];

tform3 = maketform('affine', T3);

TARDBP = tformfwd(preTARDBP, tform3)
VDAC = tformfwd(preVDAC, tform3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create mask based on oligo-dT488 image

% Read oligodT .tif file
im_mask = imread('___.tif');

% Sets threshold to only consider pixels with intensity over 2900 (depends on image), creates logical array 
BW = (im_mask > 2500); 

% Despeckle
BW2 = bwareaopen(BW, 5000);

% Find size of binary image = size of cell
size = sum(BW2(:)) 

% Get boundary / perimeter of mask and its X Y coordinates
structBoundaries = bwboundaries(BW2);

% Get n by 2 array of x,y coordinates
premaskxy=structBoundaries{1}; 
premaskx = premaskxy(:, 2); % Columns
premasky = premaskxy(:, 1) % Rows

% Coordinates of mask outline
premaskxy = [premaskx, premasky] 

%Align mask with reference
  
ts = 7.0720;
tv = 0.9991;
tw = 0.0004;
tx = 13.6656;
ty = -0.0006;
tz = 0.9985;

T4 = [tv ty 0;
      tw tz 0;
      ts tx 1]

tform4 = maketform('affine', T4);

maskxy = tformfwd(premaskxy, tform4)
maskx = maskxy(:,1);
masky = maskxy(:,2)



% Open DAPI image (of reference) and get boundary information

im_dapi = imread('dapi_hyb0_.tif');
BWdapi = (im_dapi > 3000);
BWdapi3 = imfill(BWdapi,'holes');
BWdapi2 = bwareaopen(BWdapi3, 9000),% Despeckle
dapiboundaries = bwboundaries(BWdapi2);
dxy=dapiboundaries{1}; % Get n by 2 array of x,y coordinates
dapix = dxy(:, 2); % Columns
dapiy = dxy(:, 1) % Rows
dapixy = [dapix, dapiy]


figure;
plot(maskx,masky); % mask
axis equal;

hold on;
plot(dapix,dapiy,'green');
saveas(gcf,'boundaries.png');

answer = questdlg('Close?')
switch answer
case 'Yes'
close;
end

% Find all coordinates within mask and subset for each gene

%ANTXR
Ax = A.x
Ay = A.y
inA = inpolygon(Ax, Ay, maskx, masky) 
Ainx = Ax(inA)
Ainy = Ay(inA)
Ainpoints = [Ainx, Ainy]

%NUP155
Nx = NUP155(:,1);
Ny = NUP155(:,2);
inN = inpolygon(Nx, Ny, maskx, masky) 
Ninx = Nx(inN)
Niny = Ny(inN)
Ninpoints = [Ninx, Niny]

%MYC
Mx = MYC(:,1);
My = MYC(:,2)
inM = inpolygon(Mx, My, maskx, masky) 
Minx = Mx(inM)
Miny = My(inM)
Minpoints = [Minx, Miny]

%PGRMC
Px = PGRMC(:,1);
Py = PGRMC(:,2)
inP = inpolygon(Px, Py, maskx, masky) 
Pinx = Px(inP)
Piny = Py(inP)
Pinpoints = [Pinx, Piny]

%RAVER
RVx = RAVER(:,1)
RVy = RAVER(:,2)
inRV = inpolygon(RVx, RVy, maskx, masky) 
RVinx = RVx(inRV)
RViny = RVy(inRV)
RVinpoints = [RVinx, RViny]

%RBFOX
RBx = RBFOX(:,1)
RBy = RBFOX(:,2)
inRB = inpolygon(RBx, RBy, maskx, masky) 
RBinx = RBx(inRB)
RBiny = RBy(inRB)
RBinpoints = [RBinx, RBiny]

%TARDBP
Tx = TARDBP(:,1)
Ty = TARDBP(:,2)
inT = inpolygon(Tx, Ty, maskx, masky) 
Tinx = Tx(inT)
Tiny = Ty(inT)
Tinpoints = [Tinx, Tiny]

%VDAC
Vx = VDAC(:,1)
Vy = VDAC(:,2)
inV = inpolygon(Vx, Vy, maskx, masky) 
Vinx = Vx(inV)
Viny = Vy(inV)
Vinpoints = [Vinx, Viny]

%Plot and save
figure
plot(maskx,masky,'Color',[0.4940 0.1840 0.5560]) % mask
axis equal
hold on
plot(dapix,dapiy,'Color', [0.6350 0.0780 0.1840])
hold on
plot(Ainx,Ainy,'b.')
hold on
plot(Ninx,Niny,'k.')
hold on
plot(Minx,Miny,'r.')
hold on
plot(Pinx,Piny,'g.')
hold on
plot(RBinx,RBiny,'y.')
hold on
plot(RVinx,RViny,'c.')
hold on
plot(Tinx,Tiny,'m.')
hold on
plot(Vinx,Viny,'.','Color', [0.8500, 0.3250, 0.0980])
savefig('masked_with_genes.fig');
saveas(gcf,'masked_with_genes.png');

answer = questdlg('Close?')
switch answer
case 'Yes'
close;
end

% ANTXR
[idx.am,d.am] = knnsearch(maskxy, Ainpoints, 'k', 1)
[idx.ad,d.ad] = knnsearch(dapixy, Ainpoints, 'k', 1)
AntMem = d.am
AntNuc = d.ad

%NUP
[idx.nm,d.nm] = knnsearch(maskxy, Ninpoints, 'k', 1)
[idx.nd,d.nd] = knnsearch(dapixy, Ninpoints, 'k', 1)
NupMem = d.nm
NupNuc = d.nd

%MYC
[idx.mm,d.mm] = knnsearch(maskxy, Minpoints, 'k', 1)
[idx.md,d.md] = knnsearch(dapixy, Minpoints, 'k', 1)
MycMem = d.mm
MycNuc = d.md

%PG
[idx.pm,d.pm] = knnsearch(maskxy, Pinpoints, 'k', 1)
[idx.pd,d.pd] = knnsearch(dapixy, Pinpoints, 'k', 1)
PgMem = d.pm
PgNuc = d.pd

%RAVER
[idx.rvm,d.rvm] = knnsearch(maskxy, RVinpoints, 'k', 1)
[idx.rvd,d.rvd] = knnsearch(dapixy, RVinpoints, 'k', 1)
RVMem = d.rvm
RVNuc = d.rvd

%RBFOX
[idx.rbm,d.rbm] = knnsearch(maskxy, RBinpoints, 'k', 1)
[idx.rbd,d.rbd] = knnsearch(dapixy, RBinpoints, 'k', 1)
RBMem = d.rbm
RBNuc = d.rbd

%TARDBP
[idx.tm,d.tm] = knnsearch(maskxy, Tinpoints, 'k', 1)
[idx.td,d.td] = knnsearch(dapixy, Tinpoints, 'k', 1)
TbpMem = d.tm
TbpNuc = d.td

%VDAC
[idx.vm,d.vm] = knnsearch(maskxy, Vinpoints, 'k', 1)
[idx.vd,d.vd] = knnsearch(dapixy, Vinpoints, 'k', 1)
VdcMem = d.vm
VdcNuc = d.vd

membrane.ant = mean(d.am(:));
membrane.nup = mean(d.nm(:));
membrane.myc = mean(d.mm(:));
membrane.pg = mean(d.pm(:));
membrane.rv = mean(d.rvm(:));
membrane.rb = mean(d.rbm(:));
membrane.tbp = mean(d.tm(:));
membrane.vdc = mean(d.vm(:));

writetable(struct2table(membrane), 'mean_distances.xlsx', 'Sheet',1);


dapi.ant = mean(d.ad(:));
dapi.nup = mean(d.nd(:));
dapi.myc = mean(d.md(:));
dapi.pg = mean(d.pd(:));
dapi.rv = mean(d.rvd(:));
dapi.rb = mean(d.rbd(:));
dapi.tbp = mean(d.td(:));
dapi.vdc = mean(d.vd(:));

writetable(struct2table(dapi), 'mean_distances.xlsx', 'Sheet',2);



count.cellsize = sum(BW2(:) == 1); % Number of 1s in the binary mask image => size of cell
count.ANT = numel(Ainpoints)/2;
count.MYC = numel(Minpoints)/2;
count.NUP = numel(Ninpoints)/2;
count.PGC = numel(Pinpoints)/2;
count.RB = numel(RBinpoints)/2;
count.RV = numel(RVinpoints)/2;
count.TBP = numel(Tinpoints)/2;
count.VDC = numel(Vinpoints)/2;

writetable(struct2table(count), 'spotcount.xlsx', 'Sheet',1);


% Distances between cell boundary and nuclear boundary 
filename = '8genes_in_points.xlsx';
Ainpoints = [Ainx, Ainy];
tA = table(Ainpoints);
writetable(tA,filename,'Sheet',1);

Minpoints = [Minx, Miny];
tM = table(Minpoints);
writetable(tM,filename,'Sheet',2);

Ninpoints = [Ninx, Niny];
tN = table(Ninpoints);
writetable(tN,filename,'Sheet',3);


Pinpoints = [Pinx, Piny];
tP = table(Pinpoints);
writetable(tP,filename,'Sheet',4);

RBinpoints = [RBinx, RBiny];
tRB = table(RBinpoints);
writetable(tRB,filename,'Sheet',5);


RVinpoints = [RVinx, RViny];
tRV = table(RVinpoints);
writetable(tRV,filename,'Sheet',6);


Tinpoints = [Tinx, Tiny];
tT = table(Tinpoints);
writetable(tT,filename,'Sheet',7);


Vinpoints = [Vinx, Viny];
tV = table(Vinpoints);
writetable(tV,filename,'Sheet',8);


dapixy = [dapix, dapiy]
tD = table(dapixy);
writetable(tD,filename,'Sheet',9);


maskxy = [maskx, masky]
tMASK = table(maskxy);
writetable(tMASK,filename,'Sheet',10);


%%
filenamed = 'bound_distances.xlsx';
% ANTXR
AntMem = d.am
AntNuc = d.ad
tAD = table(AntMem,AntNuc)
writetable(tAD,filenamed,'Sheet',1);

%MYC
MycMem = d.mm
MycNuc = d.md
tMD = table(MycMem,MycNuc)
writetable(tMD,filenamed,'Sheet',2);

%NUP
NupMem = d.nm
NupNuc = d.nd
tND = table(NupMem,NupNuc)
writetable(tND,filenamed,'Sheet',3);

%PG
PgMem = d.pm
PgNuc = d.pd
tPD = table(PgMem,PgNuc)
writetable(tPD,filenamed,'Sheet',4);

%RBFOX
RBMem = d.rbm
RBNuc = d.rbd
tRBD = table(RBMem,RBNuc)
writetable(tRBD,filenamed,'Sheet',5);

%RAVER
RVMem = d.rvm
RVNuc = d.rvd
tRVD = table(RVMem,RVNuc)
writetable(tRVD,filenamed,'Sheet',6);

%TARDBP
TbpMem = d.tm
TbpNuc = d.td
tTD = table(TbpMem,TbpNuc)
writetable(tTD,filenamed,'Sheet',7);

%VDAC
VdcMem = d.vm
VdcNuc = d.vd
tVD = table(VdcMem,VdcNuc)
writetable(tVD,filenamed,'Sheet',8);


%%
filenameNORM = 'normalised_distances.xlsx';

normAntMem = (AntMem/count.cellsize)*norm
normMycMem = (MycMem/count.cellsize)*norm
normNupMem = (NupMem/count.cellsize)*norm
normPgMem = (PgMem/count.cellsize)*norm
normRBMem = (RBMem/count.cellsize)*norm
normRVMem = (RVMem/count.cellsize)*norm
normTbpMem = (TbpMem/count.cellsize)*norm
normVdcMem = (VdcMem/count.cellsize)*norm




normAntNuc = (AntNuc/count.cellsize)*norm
normMycNuc = (MycNuc/count.cellsize)*norm
normNupNuc = (NupNuc/count.cellsize)*norm
normPgNuc = (PgNuc/count.cellsize)*norm
normRBNuc = (RBNuc/count.cellsize)*norm
normRVNuc = (RVNuc/count.cellsize)*norm
normTbpNuc = (TbpNuc/count.cellsize)*norm
normVdcNuc = (VdcNuc/count.cellsize)*norm




fullPath=fullfile(['/.../', 'analysismem.xlsx'])

tantmem = table(normAntMem);
writetable(tantmem,fullPath,'Sheet',1,'Range','A1');
writetable(tantmem,fullPath,'Sheet',2,'Range','A1','WriteVariableNames',false);
tmycmem = table(normMycMem);
writetable(tmycmem,fullPath,'Sheet',3,'Range','A1');
writetable(tmycmem,fullPath,'Sheet',4,'Range','A1','WriteVariableNames',false);
tnupmem = table(normNupMem);
writetable(tnupmem,fullPath,'Sheet',5,'Range','A1');
writetable(tnupmem,fullPath,'Sheet',6,'Range','A1','WriteVariableNames',false);
tpgmem = table(normPgMem);
writetable(tpgmem,fullPath,'Sheet',7,'Range','A1');
writetable(tpgmem,fullPath,'Sheet',8,'Range','A1','WriteVariableNames',false);
trvmem = table(normRVMem);
writetable(trvmem,fullPath,'Sheet',9,'Range','A1');
writetable(trvmem,fullPath,'Sheet',10,'Range','A1','WriteVariableNames',false);
trbmem = table(normRBMem);
writetable(trbmem,fullPath,'Sheet',11,'Range','A1');
writetable(trbmem,fullPath,'Sheet',12,'Range','A1','WriteVariableNames',false);
ttbpmem = table(normTbpMem);
writetable(ttbpmem,fullPath,'Sheet',13,'Range','A1');
writetable(ttbpmem,fullPath,'Sheet',14,'Range','A1','WriteVariableNames',false);
tvdcmem = table(normVdcMem);
writetable(tvdcmem,fullPath,'Sheet',15,'Range','A1');
writetable(tvdcmem,fullPath,'Sheet',16,'Range','A1','WriteVariableNames',false);




fullPathnuc = fullfile(['/.../', 'analysisnuc.xlsx'])

tantnuc = table(normAntNuc);
writetable(tantnuc,fullPathnuc,'Sheet',1,'Range','A1');
writetable(tantnuc,fullPathnuc,'Sheet',2,'Range','A1','WriteVariableNames',false);
tmycnuc = table(normMycNuc);
writetable(tmycnuc,fullPathnuc,'Sheet',3,'Range','A1');
writetable(tmycnuc,fullPathnuc,'Sheet',4,'Range','A1','WriteVariableNames',false);
tnupnuc = table(normNupNuc);
writetable(tnupnuc,fullPathnuc,'Sheet',5,'Range','A1');
writetable(tnupnuc,fullPathnuc,'Sheet',6,'Range','A1','WriteVariableNames',false);
tpgnuc = table(normPgNuc);
writetable(tpgnuc,fullPathnuc,'Sheet',7,'Range','A1');
writetable(tpgnuc,fullPathnuc,'Sheet',8,'Range','A1','WriteVariableNames',false);
trvnuc = table(normRVNuc);
writetable(trvnuc,fullPathnuc,'Sheet',9,'Range','A1');
writetable(trvnuc,fullPathnuc,'Sheet',10,'Range','A1','WriteVariableNames',false);
trbnuc = table(normRBNuc);
writetable(trbnuc,fullPathnuc,'Sheet',11,'Range','A1');
writetable(trbnuc,fullPathnuc,'Sheet',12,'Range','A1','WriteVariableNames',false);
ttbpnuc = table(normTbpNuc);
writetable(ttbpnuc,fullPathnuc,'Sheet',13,'Range','A1');
writetable(ttbpnuc,fullPathnuc,'Sheet',14,'Range','A1','WriteVariableNames',false);
tvdcnuc = table(normVdcNuc);
writetable(tvdcnuc,fullPathnuc,'Sheet',15,'Range','A1');
writetable(tvdcnuc,fullPathnuc,'Sheet',16,'Range','A1','WriteVariableNames',false);


