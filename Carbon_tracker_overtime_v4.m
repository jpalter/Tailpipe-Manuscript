%Calculate_fluxes_forTable2_Figure1and3.m

%A newer version of the code made to eliminate some bug introduced in v2
%while responding to reviewers

%the first part of this code loads everything in and is very time
%consuming.  I saved output at the end of this if-loop

if iteration == 0

    %%%% Loading data, defining directories dir 
    % First, load a single CarbonTracker file to get size of arrays, and
    % stable orography file It has one day of data at 8 time steps
    this_CT = netcdf_load('CT2019B.molefrac_glb3x2_2000-01-11.nc');
    lon = this_CT.longitude;
    lat = this_CT.latitude;

    %Use the orography variable to find land, later used to mask it
    oro = this_CT.orography;
    orot = repmat(oro,[1 1 365]);
    ind_land = orot > 10;
    [lati, loni] = meshgrid(lat,lon);
    cmphr2mpyr = 24*365/100;

    %Define directory where surface-only CT data are stored
    data_dir =  '~/Downloads/CT2019B/'
    %data_dir='/Users/jaimepalter/Library/CloudStorage/GoogleDrive-jpalter@uri.edu/Shared drives/Palter Lab/Data/Carbon Tracker/CT2019B/CT_surface_only';
    dire = dir(data_dir);
    dire = dire(3:end);

    %load and prepare Landschutzer data:
    %original manuscript used v2020
    %lan = netcdf_load('../../Lanschutzer SOM-FFN/SOMFFNv2020/spco2_MPI-SOM_FFN_v2020.nc');
    %after revisions, changed to v2022, the first one to have monthly
    %gridded product with Arctic, some marginal seas, and with consistency
    %all the way to coastlines
    lan = netcdf_load('../../Lanschutzer SOM-FFN/SOMFFNv2022/MPI_SOM_FFN_2022_NCEI_OCADS.nc');
    sol = lan.sol;
    sol(sol == max(sol(:))) = NaN; 
    kw_lan = lan.kw;%solubility [mol/m3uatm]
    spco2_smoothed = lan.spco2_smoothed; %smoothed surface ocean pCO2 [uatm]
    spco2_smoothed(spco2_smoothed==max(spco2_smoothed(:))) = NaN;
    ltime = double(lan.time);%time in seconds since 2000-01-01-00:00 (monthly resolution)
    %ldate = lan.date; year = ldate(1,:); month = ldate(2,:); %from v2020
    ltime_indecyears = 2000 + ltime/(60*60*24*365.25);
    year = floor(ltime_indecyears);
    years = 2000:1:2018;
    lonl = lan.lon;
    latl = lan.lat;
    [latil, lonil]= meshgrid(latl,lonl);
    latil = double(latil);
    lonil = double(lonil);

    %% Loop over all files to make initial calculations
    %Initialize arrays
    F_full = NaN(360,180,365,length(years));
    F_xco2zm = NaN(360,180,365,length(years));
    F_full_monthmean = NaN(360,180,365,length(years));
    pCO2G_all = NaN(120,90,365*19);
    pCO2_zmxco2_all =  NaN(120,90,365*19);
    xco2_noland_all = NaN(120,90,365*19);
    u_noland_all = NaN(120,90,365*19);
    v_noland_all = NaN(120,90,365*19);
    Patm_all = NaN(120,90,365*19);
    PH2O_all = NaN(120,90,365*19);
    %wind_speed_noland = NaN(120,90,365*19);
    kw_CT_all = NaN(120,90,365*19);
    count = 1;

     for ii =1:length(dire)  %loop over all files
        tic
        this_file =  dire(ii).name;
        %load all CT variables (surface co2, T,P, u,v specific_humidity -
        %includes leap days!):
        load(strcat(data_dir,'/',this_file)); 
        %calculate pCO2 = xCO2 * (Patm - pH2O) where Patm is
        %SLP in atmospheres and pH2O is partial pressure of water vapor in air
        %Calculate water vapor pressure with
        %Buck_pressure = 0.61121 * e^[(18.678 - (T / 234.5)) * (T / (257.14 + T))]
        %where T is expressed in °C and P in kPa.
        T(ind_land) = NaN;   %Surface Air Temp [K] 
        T = T -  272.15;

        P(ind_land) = NaN;
        specific_humidity(ind_land) = NaN;

        u(ind_land) = NaN;
        v(ind_land) = NaN;
        co2(ind_land) = NaN;
        xco2_noland = co2;
        [m,n,k] = size(xco2_noland);

        xco2_noland_all(1:m,1:n,count:count+k-1) = xco2_noland;
        u_noland_all(1:m,1:n,count:count+k-1) = u; %units m/s
        v_noland_all(1:m,1:n,count:count+k-1) = v;%units m/s
        
        %wind_speed_noland(1:m,1:n,count:count+k-1) = sqrt(u.^2 + v.^2);

        %Simplified Kw*solubility*constants from Wanninkof 2014 - units
        %will work out to mol/m2year when kw is multiplied by delta(pCO2).
        %The simplification increases uncertainty, but is necessary in the
        %absence of SST, which is not provided either by CarbonTracker or
        %Landschutzer (SAT is not a good enough approximation because of
        %very cold air temps relative to SST at high latitudes)
        kw_CT_all(1:m,1:n,count:count+k-1) = 7.7e-4*(u.^2 +v.^2);
        kw_CT_year = 7.7e-4*(u.^2 +v.^2);

         Buck_H2Opressure = 0.61121 * exp((18.678 - (T / 234.5)) .* ...
            (T ./ (257.14 + T)));
        c = P / 1.01325e5;  %convert atm press in Pa to atm
        PH2O = Buck_H2Opressure *1000 / 1.01325e5; %convert kPa to atm
        Patm = P / 1.01325e5;  %convert atm press in Pa to atm
        Patm_all(1:m,1:n,count:count+k-1) = Patm;
        PH2O_all(1:m,1:n,count:count+k-1) = PH2O;
        %full pCO2 field for year(ii) size(120x90x365)
        pCO2G = co2 .* (Patm - PH2O);
        [m,n,k] = size(pCO2G);
        pCO2G_all(1:m,1:n,count:count+k-1) = pCO2G;

        %Now calculate pCO2 with zonal mean xCO2 - note it's the zonal mean
        %only over the ocean, to be analogous with the NOAA Marine Boundary Layer
        %Zonal mean product
        zm_xco2_noland = squeeze(nanmean(xco2_noland,1));
        zm_xco2_noland = permute(zm_xco2_noland,[3 1 2]);
        xco2_zm = repmat(zm_xco2_noland,120,1,1);
        Patm_noland = Patm;
        Patm_noland(orot>10)=NaN;

        %pCO2 calculated with zonal mean xCO2 for year (ii) size(120x90x365)
        pCO2_zmxco2 = xco2_zm .* (Patm_noland - PH2O);
    
        month = ceil(12*(ltime_indecyears - year));
        pCO2_zmxco2_all(1:m,1:n,count:count+k-1) = pCO2_zmxco2;

        count = count + k;  %get all days in year (365 or 366) from 3rd dimension

        %now, find Landschutzer surface ocean pCO2 in same year to calculate
        %fluxes:
        indyear = year == years(ii); 
    
        spco2_year= spco2_smoothed(:,:,indyear);

        %replicate monthly kw, sol, spco2 for each day of the month (deal
        %with leap years)
        if years(ii) == 2000
            days_s = daysinmonths('leap'); %2000 was not a leap year!
        elseif mod(years(ii),4) ~= 0 %nonleap years
            days_s = daysinmonths('leap'); %the function is backwards "leap" gives 28 days in feb "noleap" 29!
        elseif mod(years(ii),4) == 0
            days_s = daysinmonths('noleap');
        end

        end_day = cumsum(days_s);
        start_day = [1 (end_day(1:end-1)+1)];
        spco2_dayyear = NaN(360,180,sum(days_s)); %initalize array
        %month-mean kw used to check if sub-monthly variability in winds might rectify
        %into fluxes
        kw_CT_monthmean= NaN(120,90,sum(days_s)); 
        for kk = 1:12
            kw_CT_monthmean(:,:,start_day(kk):end_day(kk)) = ...
               repmat(mean(kw_CT_year(:,:,start_day(kk):end_day(kk)),3),...
               [1 1 days_s(kk)]);  
               spco2_dayyear(:,:,start_day(kk):end_day(kk)) = repmat(spco2_year(:,:,kk),[1 1 days_s(kk)]);
        end
       
        %interpolate atmospheric pCO2 AND kw_CT to ocean grid
        pCO2Gi = NaN(size(spco2_dayyear));
        pCO2_zmxco2i = NaN(size(spco2_dayyear));
        kw_CTi = NaN(size(spco2_dayyear));
        kw_CT_monthmeani= NaN(size(spco2_dayyear));
         for jj = 1:k  %k=365 or 366, depending on the year
            this_pc = double(pCO2G(:,:,jj));
            this_pczm = double(pCO2_zmxco2(:,:,jj));
            this_kwCT_monthmean = double(kw_CT_monthmean(:,:,jj));%double(kw_CT_year(:,:,jj));
            this_kwCT = double(kw_CT_year(:,:,jj));
            %this_wind_speed = wind_speed_noland(:,:,jj);
            pCO2Gi(:,:,jj) = griddata(loni(~isnan(this_pc)),lati(~isnan(this_pc)),...
                this_pc(~isnan(this_pc)),lonil,latil);
            pCO2_zmxco2i(:,:,jj) = griddata(loni(~isnan(this_pczm)),...
                lati(~isnan(this_pczm)),this_pczm(~isnan(this_pczm)),lonil,latil);
            kw_CTi(:,:,jj) = griddata(loni(~isnan(this_pc)),...
                lati(~isnan(this_pc)),this_kwCT(~isnan(this_pc)),lonil,latil);
            kw_CT_monthmeani(:,:,jj) = griddata(loni(~isnan(this_pc)),...
                lati(~isnan(this_pc)),this_kwCT_monthmean(~isnan(this_pc)),lonil,latil);
            %wind_speedi = griddata(loni(~isnan(this_pczm)),...
            %    lati(~isnan(this_pczm)),this_wind_speed(~isnan(this_pczm)),lonil,latil);
         end %end looping through all days in the year
       
        %calculate fluxes mol m-2 year-1 (these are size 360x180x365xlength(years))
        F_full(1:360,1:180,1:sum(days_s),ii) = kw_CTi.*...
            (spco2_dayyear - pCO2Gi);
        F_xco2zm(1:360,1:180,1:sum(days_s),ii) = kw_CTi.*...
            (spco2_dayyear - pCO2_zmxco2i);
        F_full_monthmean(1:360,1:180,1:sum(days_s),ii) = kw_CT_monthmeani.*...
            (spco2_dayyear - pCO2Gi);
        ii
        toc
     end  %end looping through all CT directories

end  %end iteration signifying first iteration through code

%% Plotting and other comparisons of "full flux" with "calculated with zonal mean xCO2"
% Now also comparing CarbonTracker kw (from daily winds) with Landschutzer
if iteration == 2
%load  ~/'OneDrive - University of Rhode Island'/Tailpipe/Carbon_tracker_fluxes.mat
%load F_full.mat  %the saved output already converted using cmphr2mpyr, so

 F_fullpm =nanmean(F_full,4);
 F_fullpm = nanmean(F_fullpm,3);
 Fzm_pm = nanmean(F_xco2zm,4);
 Fzm_pm = nanmean(Fzm_pm,3);


 %%%To close the gap at the prime meridion:
latil2 = [ latil(end,:);latil;latil(end,:) ];
lonil2 = [lonil(1,:)-0.5; lonil;lonil(end,:)+0.5];
F_fullpm2 = F_fullpm;
F_fullpm2(1,:)=F_fullpm(2,:);  %NaNs in first and final position need to be replaced
F_fullpm2(end,:) = F_fullpm(end-1,:);
F_fullpm2 = [F_fullpm2(1,:); F_fullpm2; F_fullpm2(end,:)];

%%%
figure('Position', [100, 200, 800, 600])
%axesm('MapProjection','mollweid','MapLonLimit',[0 360]);framem; gridm

 axesm('robinson', 'frame', 'on','origin',[0 180 0]);
 gridm
%ax = worldmap('world');
%setm(ax,"Origin",[0 180])
title('a) Average air-sea CO_2 flux (mol m^-^2 year^-^1)')
%contourfm(latil2,lonil2,F_fullpm2,'linestyle','none');%shading flat;colorbar%shading flat;
pcolorm(double(latil2),double(lonil2),F_fullpm2);shading flat;colorbar;
cmocean('balance',20);
colorbar
hold on
%[c h] = contourm(latiG,loniG,std(surf_co2_2019missionG,1,3),[0.5:0.5:2.5],'w')
caxis([-6 6])
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeColor','k')%[.7 .7 .7])

%Other panels haven't changed from v1 of carbon_tracker_overtime (moved and
%commented for completeness):
%select ocean pco2 2000-2018 and average
spco2_gy = spco2_smoothed(:,:,217:end);
spco2_tm = mean(spco2_gy,3);
spco2_std = std(spco2_gy,1,3);
figure('Position', [100, 200, 800, 600])
spco2_tm2 = [spco2_tm(1,:); spco2_tm; spco2_tm(end,:)];
ax = worldmap('world');
setm(ax,"Origin",[0 180 0])
%axesm('MapProjection','eqdcylin','MapLatLimit',[10 50],'MapLonLimit',[-85 10])
title('b) Ocean surface pCO_2 (\muatm)')
pcolorm(latil2,lonil2,spco2_tm2);shading flat;colorbar%shading flat;
hold on;
contourm(latil,lonil,double(spco2_std),[20 20],'k');
cmocean('haline',20);
colorbar
hold on
%[c h] = contourm(latiG,loniG,std(surf_co2_2019missionG,1,3),[0.5:0.5:2.5],'w')
caxis([310 480])
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeColor','k')%[.7 .7 .7])

%rcepeat for atmospheric pCO2
pCO2G_allmean = nanmean(pCO2G_all,3);%NaNs in final year
%monthly standard deviation requires first averaging days into months
dd = daysinmonths('leap');
dd = cumsum(dd);
dd = repmat(dd,20,1);
dd = dd(:);
dd = [1; dd];
for ii = 1:length(dd)-1
    this_month = pCO2G_all(:,:,dd(ii):dd(ii+1));
    pCO2G_allmonths(1:120,1:90,ii) = mean(this_month,3);
end
%now take the monthly standard deviation
pCO2G_mstd = nanstd(pCO2G_allmonths,0,3);
%
lati2 = [lati; lati(end,:)];
loni2 = [loni;loni(1,:)+0.5];

figure('Position', [100, 200, 800, 600])
ax = worldmap('world');
setm(ax,"Origin",[0 180 0])
%axesm('MapProjection','eqdcylin','MapLatLimit',[10 50],'MapLonLimit',[-85 10])
title('c) \DeltapCO_2 (\muatm)')
%put atmospheric pCO2 on landschutzer grid
pCO2G_allmean2 = griddata(loni(:),lati(:),pCO2G_allmean(:),lonil,latil);
%close gap at prime meridion
pCO2G_allmean2 = [pCO2G_allmean2(end-2:end-1,:); pCO2G_allmean2(2:end-1,:); pCO2G_allmean2(end-2:end-1,:)];
pcolorm(lati2,loni2,spco2_tm2-pCO2G_allmean2);shading flat;colorbar%shading flat;
hold on;
cmocean('balance',20);
caxis([-100 100])
colorbar
hold on
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5],'EdgeColor','k')%[.7 .7 .7])

figure('Position', [100, 200, 800, 600])
ax = worldmap('world');
setm(ax,"Origin",[0 180 0])
%axesm('MapProjection','eqdcylin','MapLatLimit',[10 50],'MapLonLimit',[-85 10])
title('c) Atmospheric pCO_2 (\muatm)')
%close gap at prime meridion

pCO2G_allmean2 = [pCO2G_allmean; pCO2G_allmean(end,:)];
pcolorm(lati2,loni2,pCO2G_allmean2);shading flat;colorbar%shading flat;
hold on;
%contourm(latil,lonil,double(spco2_std),[20 20],'k');
cmocean('haline',20);
colorbar
hold on
contourm(lati,loni,pCO2G_mstd,[10 10],'k');
caxis([310 480]);
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5],'EdgeColor','k')%[.7 .7 .7])



%repeat for atmospheric xCO2 anomaly from zonal mean
%first remove empty months at the end
xco2_noland_all(:,:,6940:end) = [];
u_noland_all(:,:,6940:end) = [];
%wind_speed_noland(:,:,6940:end)=[];
zm_xco2_all = squeeze(nanmean(xco2_noland_all,1));
zm_xco2_all = permute(zm_xco2_all,[3 1 2]);
xco2_zm_all = repmat(zm_xco2_all,120,1,1);
xco2_anom = xco2_noland_all - xco2_zm_all;
dd = daysinmonths('leap');
dd = cumsum(dd);
dd = repmat(dd,20,1);
dd = dd(:);
dd = [1; dd];
%monthly standard deviation requires first averaging days into months
for ii = 1:length(dd)-1
    this_month = xco2_anom(:,:,dd(ii):dd(ii+1));
    xCO2_allmonths(1:120,1:90,ii) = mean(this_month,3);
end
%now take the monthly standard deviation
xCO2anom_mstd = nanstd(xCO2_allmonths,0,3);

figure('Position', [100, 200, 800, 600])
ax = worldmap('world');
setm(ax,"Origin",[0 180 0])
%axesm('MapProjection','eqdcylin','MapLatLimit',[10 50],'MapLonLimit',[-85 10])
title('d) Atmospheric xCO_2 anomaly from zonal mean (\muatm)')
xco2_anom_mean = mean(xco2_anom,3);
xco2_anom_mean2 = [xco2_anom_mean; xco2_anom_mean(end,:)];
contourfm(lati2,loni2,xco2_anom_mean2,[-4:0.5:4],'linestyle','none');%shading flat;colorbar%shading flat;
hold on;
%contourm(latil,lonil,double(spco2_std),[20 20],'k');
cmocean('balance',20);
colorbar
hold on
contourm(lati,loni,xCO2anom_mstd,[2 2],'k');
caxis([-3.5 3.5])
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5],'EdgeColor','k')%[.7 .7 .7])


%NORTHERN HEMISPHERE flux plots
tailpipe_contour = -0.05;
figure;
subplot(3,1,1);
ax = axesm('MapProjection','mercator','MapLatLimit',[5 55],'MapLonLimit',[-180 180]);
setm(ax,"Origin",[0 180 0])

land = readgeotable("landareas.shp");
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.25,'EdgeColor',[0.5 0.5 0.5]);%[.7 .7 .7])
pcolorm(latil2,lonil2,F_fullpm2);shading flat;colorbar
tightmap
caxis([-6 6])
geoshow('landareas.shp','FaceColor','k','FaceAlpha',0.25,'EdgeColor',[0.5 0.5 0.5])%[.7 .7 .7])
cmocean('balance',20);
colorbar
title('Mean flux 2000-2018 (mol m^-^2 year^-^1)')

%close gap at prime meridion
Fzm_pm2 = Fzm_pm;
Fzm_pm2(1,:) = Fzm_pm(2,:);
Fzm_pm2(end,:) = Fzm_pm(end,:);
Fzm_pm2 = [Fzm_pm2(1,:); Fzm_pm2;Fzm_pm2(end,:)];

subplot(3,1,2);
ax = axesm('MapProjection','mercator','MapLatLimit',[5 55],'MapLonLimit',[-180 180]);
setm(ax,"Origin",[0 180 0])
land = readgeotable("landareas.shp");
pcolorm(latil2,lonil2,F_fullpm2 - Fzm_pm2);shading flat;colorbar
hold on
perc_diff = 100*(F_fullpm-Fzm_pm)./abs(F_fullpm);
indhigh = abs(perc_diff) >= 10 & abs(F_fullpm) >=0.5;
hold on;
%plotm(latil(indhigh),lonil(indhigh),'k+','markersize',0.25);%,'color',[0.5 0.5 0.5])
contourm(latil,lonil,F_fullpm - Fzm_pm,[tailpipe_contour tailpipe_contour],'k','linewidth',1.5);shading flat;colorbar

tightmap
geoshow('landareas.shp','FaceColor','k','FaceAlpha',0.25,'EdgeColor',[0.5 0.5 0.5])%[.7 .7 .7])
caxis([-.5 .5])
cmocean('balance',20);
colorbar
title('Mean flux difference (mol m^-^2) year^-^1')

%repeat for MAX difference!
%figure;
subplot(3,1,3);

ax = axesm('MapProjection','mercator','MapLatLimit',[5 55],'MapLonLimit',[-180 180]);
setm(ax,"Origin",[0 180 0])
land = readgeotable("landareas.shp");
F_diff_all = (F_full - F_xco2zm);
[Fdiff_alltime_min indmin] = min(F_diff_all,[],[3 4],'linear');
F_full_at_diffmin = F_full(indmin);
perc_diff_atmin = 100*Fdiff_alltime_min./F_full_at_diffmin;

%close gap at prime meridion
Fdiff_alltime_min2 = Fdiff_alltime_min;
Fdiff_alltime_min2(1,:)=Fdiff_alltime_min(2,:);
Fdiff_alltime_min2(end,:) =Fdiff_alltime_min(end,:);
Fdiff_alltime_min2 = [Fdiff_alltime_min2(1,:); Fdiff_alltime_min2;Fdiff_alltime_min2(end,:)];
pcolorm(latil2,lonil2,Fdiff_alltime_min2);shading flat;colorbar
hold on
ind3 = perc_diff_atmin >30;
plotm(latil(ind3),lonil(ind3),'k.')

%contourm(latil,lonil,100*(F_fullpm-Fzm_pm)./abs(F_fullpm),[-5 -5],'k-')
%contourm(latil,lonil,100*(F_fullpm-Fzm_pm)./abs(F_fullpm),[5 5],'k')
tightmap
geoshow('landareas.shp','FaceColor','k','FaceAlpha',0.25,'EdgeColor',[0.5 0.5 0.5])%
caxis([-3 3])
cmocean('balance',20);
colorbar
title('Max daily underestimate (mmol m^-^2 year^-^1)')
mlabel('MLabelLocation',[0:60:360],'MLabelParallel','south');

%Carioca Buoy (magenta) and Sable Island (black) locations
%plotm(44.3,-66.3,'m.','markersize',15)
%plotm(43.9,-59.9,'k.','markersize',15)


%% Integrate over different areas to see Full versus ZM fluxes
%Define a "tailpipe" contour (i.e. where using zonal mean atm xCO2 would
%lead to an underestimate of ocean carbon uptake east of Asia or North
%America) and produce time series of uptake with full atm xCO2 and zonal
%mean atm xCO2

tailpipe_contour = -0.05;
ind_asia_tailpipe = F_fullpm-Fzm_pm <= tailpipe_contour & lonil <= 180 & lonil > 100 & latil >10;
hold on;plotm(latil(ind_asia_tailpipe),lonil(ind_asia_tailpipe),'k.')

ind_NA_tailpipe = F_fullpm-Fzm_pm <= tailpipe_contour & lonil <= -55 & lonil > -100 & latil >10;
plotm(latil(ind_NA_tailpipe),lonil(ind_NA_tailpipe),'k.')

ind_NH = latil>0;

%area factors for integrating
dy = 111e3*ones(size(latil));  %1deg latitude = 111e3 m
dx = 111e3.* cos(latl*pi/180);
dx = repmat(dx',[360 1]);
area_t = dx.*dy;
area_t_rep = repmat(area_t,[1 1 365]);
sum_Asia_tailpipe_full = NaN(365,19);
sum_NA_tailpipe_full = NaN(365,19);
sum_Asia_tailpipe_zm = NaN(365,19);
sum_NA_tailpipe_zm = NaN(365,19);
sum_global_full = NaN(365,19);
sum_global_zm = NaN(365,19);

ind_SCS = find(lonil > 132 & lonil<138 & latil>38 & latil<44);
ind_MAB = find(lonil > -75& lonil<-70 & latil>34 & latil<42);
sum_SCS_full = NaN(365,19);
sum_SCS_zm = NaN(365,19);
sum_MAB_full = NaN(365,19);
sum_MAB_zm = NaN(365,19);


area_asia_tailpipe_km2 = sum(area_t(ind_asia_tailpipe))/1e6;
area_NA_tailpipe_km2 = sum(area_t(ind_NA_tailpipe))/1e6;

area_SCS_km2 = sum(area_t(ind_SCS))/1e6;
area_MAB_km2 = sum(area_t(ind_MAB))/1e6;

count = 1;

for ii = 1:length(dire) %loop over all years
    for jj = 1:365 %loop over all days
        this_F_full = F_full(:,:,jj,ii).*area_t; %mol/year in each grid cell
        this_Fzm = F_xco2zm(:,:,jj,ii).*area_t; %mol/year in each grid cell
%have to use nansum for asia tailpipe because the new Landschutzer product
%has 3 cells near the coast of Asia that have data in all months but
%february!
        sum_Asia_tailpipe_full(jj,ii) = nansum(this_F_full(ind_asia_tailpipe));
        sum_NA_tailpipe_full(jj,ii) = sum(this_F_full(ind_NA_tailpipe)); %mol/year in region
        sum_Asia_tailpipe_zm(jj,ii) = nansum(this_Fzm(ind_asia_tailpipe));
        sum_NA_tailpipe_zm(jj,ii) = sum(this_Fzm(ind_NA_tailpipe));

        sum_global_full(jj,ii)= nansum(this_F_full(:));
        sum_global_zm(jj,ii) = nansum(this_Fzm(:));
        sum_NH_full(jj,ii)= nansum(this_F_full(ind_NH));
        sum_NH_zm(jj,ii) = nansum(this_Fzm(ind_NH));

        count = count +1;
    end
end

%variables starting with sums, are integrals in mol/year of the fluxes over
%each region.  They are mxn with m=365 days in year (neglecting leap) and
%n=19 years.

areas = ["Global";"NH";"NATailpipe";"AsiaTailpipe"];%"NH Subtropical Gyres",
mol2Pg = 12/1e15;
Full_flux = round(mol2Pg*[mean(sum_global_full(:));mean(sum_NH_full(:));...
    mean(sum_NA_tailpipe_full(:));mean(sum_Asia_tailpipe_full(:));],3);
ZM_flux = round(mol2Pg*[mean(sum_global_zm(:));mean(sum_NH_zm(:));...
    mean(sum_NA_tailpipe_zm(:));mean(sum_Asia_tailpipe_zm(:));],3);
Full_std = round(mol2Pg*[std(sum_global_full(:));std(sum_NH_full(:));...
    std(sum_NA_tailpipe_full(:));std(sum_Asia_tailpipe_full(:));],2);
ZM_std = round(mol2Pg*[std(sum_global_zm(:));std(sum_NH_zm(:));...
    std(sum_NA_tailpipe_zm(:));std(sum_Asia_tailpipe_zm(:));],2);
Full_minus_ZMperc = round(100*[(mean(sum_global_full(:))-mean(sum_global_zm(:)))/mean(sum_global_full(:));...
    (mean(sum_NH_full(:))-mean(sum_NH_zm(:)))/mean(sum_NH_full(:));...
    (mean(sum_NA_tailpipe_full(:))-mean(sum_NA_tailpipe_zm(:)))/ mean(sum_NA_tailpipe_full(:));...
    (nanmean(sum_Asia_tailpipe_full(:))-nanmean(sum_Asia_tailpipe_zm(:)))/ nanmean(sum_Asia_tailpipe_full(:));],2);


Table2 = table(areas,Full_flux,Full_std,ZM_flux,...
    ZM_std,Full_minus_ZMperc);



end

