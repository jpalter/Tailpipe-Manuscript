

clear all;%close all

%load saildrone2021_2022_licor.mat
load SD_Subset_for_Tailpipe_paper.mat %Instead of loading large Saildrone Files, load in just the subset needed to make Figure 2
load gebco_NE.mat

%% subset GEBCO 30 arcsec topography
% ilat_topo = find(lat_gebco >= 36 & lat_gebco <= 42.5);
% ilon_topo1 = find(lon_gebco >= -80 & lon_gebco <= -68);
% elev_subset = elev(ilon_topo1,ilat_topo);
% 
% lat_plot = lat_gebco(ilat_topo);
% lon_plot = lon_gebco(ilon_topo1);

%% Read in CT-NRT_v2023 North Atlantic ocean files
files = dir('*.nc');
for ii = 1:length(files)
    ctnrt(ii) = netcdf_load(files(ii).name);
end

% convert time to mtime
for ik = 1:length(files)
    ctnrt(ik).mtime = double(datenum('2000', 'yyyy') + ctnrt(ik).time);% / 864e2;
end

% grab the surface level data
for io = 1:length(files)
    if io == 1
        mtime = ctnrt(io).mtime;
        surf_xco2 = double(permute(squeeze(ctnrt(io).co2(:,:,1,:)),[2 1 3]));
        pressure = double(permute(squeeze(ctnrt(io).pressure(:,:,1,:)),[2 1 3]));
        specific_humidity = double(permute(squeeze(ctnrt(io).specific_humidity(:,:,1,:)),[2 1 3]));
    else
        this_surf_xco2 = double(permute(squeeze(ctnrt(io).co2(:,:,1,:)),[2 1 3]));
        this_surf_press = double(permute(squeeze(ctnrt(io).pressure(:,:,1,:)),[2 1 3]));
        this_surf_humid = double(permute(squeeze(ctnrt(io).specific_humidity(:,:,1,:)),[2 1 3]));
        mtime = cat(1,mtime,ctnrt(io).mtime);
        surf_xco2 = cat(3,surf_xco2,this_surf_xco2);
        pressure = cat(3,pressure,this_surf_press);
        specific_humidity = cat(3,specific_humidity,this_surf_humid);
    end
end
longitude = ctnrt(1).longitude;
latitude = ctnrt(1).latitude;

%ind_1091 = find(saildrone1091.mtime <= datenum(2021,12,18,23,59,59));
%ind_1091co2 = find(saildrone1091.co2_mtime <= datenum(2021,12,18,23,59,59));
%ind_1092co2 = find(saildrone1092.co2_mtime <= datenum(2021,12,18,23,59,59));
%ind_1090co2 = find(saildrone1090.co2_mtime <= datenum(2021,12,18,23,59,59));

surf_xco2_map = nanmean(surf_xco2,3);
surf_xco2_interp1091 = interp3(longitude,latitude,mtime,surf_xco2,saildrone1091.longitude(ind_1091),saildrone1091.latitude(ind_1091),saildrone1091.mtime(ind_1091));


%xco2_atm_1090 = saildrone1090.xco2_dry_air_mean_asvco2;
%xco2_atm_1090 = xco2_atm_1090(~isnan(xco2_atm_1090));
%xco2_atm_1091 = saildrone1091.xco2_dry_air_mean_asvco2;
%xco2_atm_1091 = xco2_atm_1091(~isnan(xco2_atm_1091));
%xco2_atm_1092 = saildrone1092.xco2_dry_air_mean_asvco2;
%xco2_atm_1092 = xco2_atm_1092(~isnan(xco2_atm_1092));
%xco2_atm_1092 = saildrone1092.xco2_dry_air_mean_asvco2;
%xco2_atm_1092 = xco2_atm_1092(~isnan(xco2_atm_1092));

%lat1090 = interp1(saildrone1090.mtime,saildrone1090.latitude,saildrone1090.co2_mtime);
%lon1090 = interp1(saildrone1090.mtime,saildrone1090.longitude,saildrone1090.co2_mtime);
%lat1091 = interp1(saildrone1091.mtime,saildrone1091.latitude,saildrone1091.co2_mtime);
%lon1091 = interp1(saildrone1091.mtime,saildrone1091.longitude,saildrone1091.co2_mtime);
%lat1092 = interp1(saildrone1092.mtime,saildrone1092.latitude,saildrone1092.co2_mtime);
%lon1092 = interp1(saildrone1092.mtime,saildrone1092.longitude,saildrone1092.co2_mtime);
%lat1092 = interp1(saildrone1092.mtime,saildrone1092.latitude,saildrone1092.co2_mtime);
%lon1092 = interp1(saildrone1092.mtime,saildrone1092.longitude,saildrone1092.co2_mtime);
latlim = [38.5 42.5];%[37.5 42.5];%[35.5 42.5]
lonlim =  [-73.5 -68.485];% [-75.5 -62]


[lati loni]= meshgrid(latitude,longitude);
figure;
ax = axesm('MapProjection','mercator','MapLatLimit',latlim,'MapLonLimit',...
    lonlim, 'MLabelLocation',[-73 -70],'PLabelLocation',[39 41],...
    'Grid','on','MLineLocation',2,'PLineVisible', 'on','MLineVisible', 'on','PLineLocation',5,...
    'fontsize',20,'MlabelParallel','south');
framem on; mlabel on; plabel on; gridm on
surfm(lati,loni,surf_xco2(:,:,17)');shading flat;
c = parula(10);colormap(c);
caxis([415 450]);y1=colorbar; ylabel(y1,'Atm xCO_2 [\mumol/mol]','fontsize',20)
contourm(lat_plot,lon_plot,elev_subset',[0 0],'k','linewidth',2)
scatterm(flipud(lat1090(ind_1090co2)),flipud(lon1090(ind_1090co2)),30,flipud(xco2_atm_1090(ind_1090co2)),'filled')
scatterm(flipud(lat1091(ind_1091co2)),flipud(lon1091(ind_1091co2)),30,flipud(xco2_atm_1091(ind_1091co2)),'filled')
scatterm(flipud(lat1092(ind_1092co2)),flipud(lon1092(ind_1092co2)),30,flipud(xco2_atm_1092(ind_1092co2)),'filled')
set(gca,'fontsize',20);tightmap
plotm(lat1090(ind_1090co2),lon1090(ind_1090co2),'k','linewidth',0.3)
plotm(lat1091(ind_1091co2),lon1091(ind_1091co2),'k','linewidth',0.3)
plotm(lat1092(ind_1092co2),lon1092(ind_1092co2),'k','linewidth',0.3)

datelabels = [datenum(2021,12,11);datenum(2021,12,12);datenum(2021,12,13); datenum(2021,12,14); datenum(2021,12,15)];
for ii= 1:5
    [minValue,closestIndex] = min(abs(saildrone1091.co2_mtime-datelabels(ii)));
    [minValue,closestIndex90] = min(abs(saildrone1090.co2_mtime-datelabels(ii)));
    dv = datevec(saildrone1091.co2_mtime(closestIndex));
    dv90 = datevec(saildrone1090.co2_mtime(closestIndex90));
    textm(lat1091(closestIndex),lon1091(closestIndex),...
        strcat(num2str(dv(2)),'.',num2str(dv(3))),...
        'color','w');
    if ii == 4
       textm(lat1090(closestIndex90),lon1090(closestIndex90),...
        strcat(num2str(dv90(2)),'.',num2str(dv90(3))),...
        'color','w');
    end
end
textm(39,-73,'SD1090');
textm(40,-69,'SD1091');
textm(40.5,-70,'SD1092')

figure;hold on
plot(saildrone1090.co2_mtime(ind_1090co2),xco2_atm_1090(ind_1090co2),'linewidth',2)
plot(saildrone1091.co2_mtime(ind_1091co2),xco2_atm_1091(ind_1091co2),'linewidth',2)
plot(saildrone1092.co2_mtime(ind_1092co2),xco2_atm_1092(ind_1092co2),'linewidth',2)
plot(saildrone1091.mtime(ind_1091),surf_xco2_interp1091,'k','linewidth',2)
set(gca,'fontsize',20);grid
ylabel('Atmospheric xCO_2 [\mumol/mol]','fontsize',20)
datetick('x','mm/dd','keepticks');
legend('SD-1090','SD-1091','SD-1092','CT-NRTv2023 interp')

