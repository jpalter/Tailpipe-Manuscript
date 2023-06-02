%carbon_tracker_case_studies.m

%compare fluxes calculated with CT versus other ad hoc methods as published
%in Rutherford et al. for Scotian Shelf (Carioca Buoy, using Sable Island
%Seasonal cycle) and Liu et al. 2022 in the East China Sea.

%Will need to add locations to Figure 2 in carbon_tracker_overtime_v2.m
% Sable Island 43.9N, 59.9 W; %out of tailpipe, source of atm CO2 info
% Carioca buoy 44.3N, 66.3 W;  In Tailpipe

%load  ~/'OneDrive - University of Rhode Island'/Tailpipe/Carbon_tracker_fluxes.mat
%load ~/'OneDrive - University of Rhode Island'/Tailpipe/pCO2_zmxco2_all

%some variables I'll need
%kw_CT 120x90x7300  gas exchange coefficient lon,lat,time
%ln  struct array with kw grid
%kw_year

%xco2_noland_all
%zm_xco2_noland
%spco2_smoothed 360,180,456 landschutzer surface ocean pCO2 lon,lat,month/year
%pCO2G_all 120x90x7300
%wind_speed_noland 120x90x7300 CT windspeed
%month and year %1x456 landschutzer month and year
%Flux equation = 7.7e-4.* wind_speedi.^2 .* (spco2_dayyear - pCO2Gi)*kwscale_factor;
kwscale_factor =  0.8403;  %found by dividing the 14C constraint by the unconstrained kw (found in carbon_tracker_overtime_v2.m)

CTtime = 1999:1/365:2019-(1/365);  %decimal time for carbon tracker
CTtime(6940:end) = [];
indg= ~isnan(pCO2G_all(38,67,:));
figure;subplot(2,1,1);
plot(CTtime,squeeze(pCO2G_all(38,67,indg)),'color',[0.5 0.5 0.5])
hold on
plot(CTtime,squeeze(pCO2_zmxco2_all(38,67,indg)),'r')
%Sable Island pCO2 empirical eq from Rutherford Supplement
CTtimeSI = CTtime-1951;
Sable_island_pco2_atm = 273.9904 + 1.9738*CTtimeSI + 5.3132*cos(2*pi*CTtimeSI)...
    + 4.0601*sin(2*pi*CTtimeSI) - 0.5814*cos(4*pi*CTtimeSI) - ...
    2.1740*sin(4*pi*CTtimeSI)-0.5009*cos(6*pi*CTtimeSI)+ 0.5904*sin(6*pi*CTtimeSI);
Sable_island_pco2_atm = Sable_island_pco2_atm';
hold on;
carioca_CTaCO2 = squeeze(pCO2G_all(38,67,:));
carioca_zm = squeeze(pCO2_zmxco2_all(38,67,:));

indnan = ~isnan(carioca_CTaCO2);
filt_CTaCO2 = filtfilt(ones(30,1)/30,1,carioca_CTaCO2(indnan));
indnan2 = ~isnan(carioca_zm);
filt_zm = filtfilt(ones(30,1)/30,1,carioca_zm(indnan2));

hold on; plot(CTtime(indnan),filt_CTaCO2,'k','LineWidth',2)
plot(CTtime(indnan2),filt_zm,'r','LineWidth',2)
plot(CTtime,Sable_island_pco2_atm,'b','linewidth',1);

axis([2015 2017 385 420])
%netcdf is in strange format, just import data from text download
%accessed here https://gaw.kishou.go.jp/ on 3/21/2023, saved as
%Sable_island_daily.mat
%sable_island_wsa = netcdf_load('co2_wsa_surface-insitu_20_9999-9999_hourly.nc')
load ~/'OneDrive - University of Rhode Island'/Tailpipe/Sable_island_daily.mat
si_yy = WSA(:,2);
si_mm = WSA(:,3);
si_day = WSA(:,4);
si_xco2 = WSA(:,14);
si_dd = si_yy + si_mm/12 -(1/12) + si_day/365 - (si_day/365);

%Sable Island 43.9N, 59.9 W
Patm_interp_si = interp1(CTtime,squeeze(Patm(41,67,indg)),si_dd);
PH2O_interp_si = interp1(CTtime,squeeze(PH2O(41,67,indg)),si_dd);
WAS_pCO2 =  si_xco2 .* (Patm_interp_si - PH2O_interp_si);
hold on
plot(si_dd,WAS_pCO2,'b+')

subplot(2,1,2)
%compare fluxes - but what ocaen pCO2?  Is data really in SOCAT?

%spco2_smoothed or spco2_dayyear needs to be selected at correct lon,lat
%and then interpolated for each day of month.
%nearest populated lon/lat is spco2_dayyear(114,133) (-66.5, 42.5)
%CT is daily, Landschutzer is monthly.  Interpolate monthly to daily
landate = ldate(1,:) + ldate(2,:)/12 - 1/24;
interp_spco2 = interp1(landate',squeeze(spco2_smoothed(114,133,:)),CTtime);

carioca_aCO2 = squeeze(pCO2G_all(38,67,:));
carioca_ws = squeeze(wind_speed_noland(38,67,:));
%DATE IS MONTH/12 - 1/12+ YEAR
%check units!  Does not look like rutherford!


CTdiff = (interp_spco2(indg)'- carioca_aCO2(indg));
SIdiff = (interp_spco2(indg)'- Sable_island_pco2_atm);
ZMdiff = (interp_spco2'- carioca_zm(indg));
%WRONG, but not terrible.  Do a better job later
%F_SI = squeeze(kw_CT(38,67,5841:5841+364)).* ...
 %  mean(sol_dayyear(114,133,:),3).*(SIdiff(5841:5841+364));%% Sol_dayyear is just a seasonal cycle

%plot(CTtime(indnan),localflux_wCT(indnan),'color',[0.5 0.5 0.5]);hold on
%filt_fluxCT = filtfilt(ones(30,1)/30,1,double(localflux_wCT(indnan)));
hold on; plot(CTtime,CTdiff(indg),'k');
plot(CTtime,SIdiff,'g')
plot(CTtime,ZMdiff(indg),'r')

filtdiffCT = filtfilt(ones(30,1)/30,1,double(CTdiff(indg)));
filtdiffSI = filtfilt(ones(30,1)/30,1,double(SIdiff(indg)));
filtdiffZM = filtfilt(ones(30,1)/30,1,double(ZMdiff));

plot(CTtime,filtdiffCT,'k','linewidth',2);
plot(CTtime,filtdiffSI,'g','linewidth',2);
plot(CTtime,filtdiffZM,'r','linewidth',2)

%plot(CTtime(indnan2),localflux_wzm(indnan2),'r');hold on
indnan2 = ~isnan(localflux_wzm);
filt_fluxzm = filtfilt(ones(30,1)/30,1,double(localflux_wzm(indnan2)));
hold on; plot(CTtime(indnan2),filt_fluxzm,'r','linewidth',1);

%plot(CTtime,localflux_wSI,'g');
filt_fluxSI= filtfilt(ones(30,1)/30,1,double(localflux_wSI(indnan)));
plot(CTtime(indnan),filt_fluxSI,'g','linewidth',1);

%axis([2015 2018 -3.75 0.5])
%subplot(2,1,2)
%plot(CTtime,filt_fluxCT(1:end-4)-filt_fluxzm)


latil(113,133)
lonil(113,133)
localflux_wCT = squeeze(F_full(113,133,:,:));
localflux_wzm = squeeze(F_xco2zm(113,133,:,:));
this_year = 2014 - 1999;
nanmean(localflux_wCT(:,this_year));
nanmean(localflux_wzm(:,this_year));

% From Liu et al., 2022
%The study region covers a 4° × 2° area from 28.5 to 32.5°N and from
% 122 to 124°E (the rectangle in Figure 1a), located in the outer
% Changjiang estuary of the ECS. Its western part is mainly affected by the
% Changjiang River discharge, where high pCO2 freshwater meets the highly
% productive ECSSW (Zhai et al., 2007; Yu et al., 2013). 
% 
% Atmospheric CO2
% concentrations were derived from the Tae-ahn Peninsula observation site
% (36.7376°N, 126.1328°E; Republic of Korea1), after correction for water
% vapor pressure to 100% humidity with in situ temperature and salinity
% data (Weiss and Price, 1980)
%East China Sea
%accessed here https://gml.noaa.gov/
%saved TAP_flask_CO2samples copydataonly - need to save as .mat

%closest point to averaging box, but NaNs
lati(101,61); %31N
loni(101,61)%121.5E
%need to go to
loni(102,61);%124.5 E

%Tae-ahn Peninsula (TAP flask station NOAA-ESRL CO2 stations)
lati(102,63)%35N
loni(102,63)%125.5N

subplot(2,1,2);
plot(CTtime,squeeze(pCO2G_all(102,63,indg)),'color',[0.5 0.5 0.5])
hold on
plot(CTtime,squeeze(pCO2_zmxco2_all(102,63,indg)),'r')

%import esrl data from TAP station, South KoreaT
%now at /OneDrive/Tailpipe/TAP_flaskCO2samples_withMeta.txt
TAP_xCO2 = TAP(:,11); %Mole fraction reported in units of micromol mol-1 (10-6 mol per mol of dry air); abbreviated as ppm (parts per million).
TAP_dd = TAP(:,2) + TAP(:,3)/12 - (1/12) + TAP(:,4)/365 - (1/365);
Patm_interp = interp1(CTtime,squeeze(Patm(102,63,indg)),TAP_dd);
PH2O_interp = interp1(CTtime,squeeze(PH2O(102,63,indg)),TAP_dd);
TAP_pCO2 =  TAP_xCO2 .* (Patm_interp - PH2O_interp);


filt_flux_ecs = filtfilt(ones(30,1)/30,1,squeeze(pCO2G_all(102,63,indg)));
filt_flux_ecszm = filtfilt(ones(30,1)/30,1,squeeze(pCO2_zmxco2_all(102,63,indg)));
plot(CTtime',filt_flux_ecs,'k','linewidth',2);
hold on; plot(CTtime,filt_flux_ecszm,'r','linewidth',2);



latil(305,122)
lonil(305,122)
localflux_wCT = squeeze(F_full(305,122,:,:));
localflux_wzm = squeeze(F_xco2zm(305,122,:,:));
this_year = 2014 - 1999;
nanmean(localflux_wCT(:,this_year))
nanmean(localflux_wzm(:,this_year))

%for insets:
figure;

ax = axesm('MapProjection','mercator','MapLatLimit',[40.5 45.5],'MapLonLimit',...
    [-75.5 -58], 'MLabelLocation',[-74 -70 -66 -62],'PLabelLocation',[36 38 40 42],...
    'Grid','on','MLineLocation',2,'PLineVisible', 'on','MLineVisible', 'on','PLineLocation',5,...
    'fontsize',20,'MlabelParallel','south');
framem on; mlabel off; plabel off; gridm on
pcolorm(latil,lonil,squeeze(nanmean(F_diff_all(:,:,:,this_year),3)));
shading flat;
colorbar
geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
set(gca,'fontsize',20);tightmap
hold on
plotm(43,-67.5,'ko','markersize',10,'markerfacecolor','k');
%plotm(43.9,-60,'b+','markersize',10,'linewidth',2)
plotm(43.9,-60,'bo','markersize',10,'linewidth',2)
cmocean('balance',20);
caxis([-.6 .6])

%for insets:
figure;
ax = axesm('MapProjection','mercator','MapLatLimit',[28 40],'MapLonLimit',...
    [114 140]);%, 'MLabelLocation',[-74 -70 -66 -62],'PLabelLocation',[36 38 40 42],...
    %'Grid','on','MLineLocation',2,'PLineVisible', 'on','MLineVisible', 'on','PLineLocation',5,...
   % 'fontsize',20,'MlabelParallel','south');
framem on; mlabel off; plabel off; gridm on
pcolorm(latil,lonil,squeeze(nanmean(F_diff_all(:,:,:,this_year),3)));
shading flat;
colorbar
geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
set(gca,'fontsize',20);tightmap
hold on
plotm(31.5,124.5,'ko','linewidth',4,'markersize',10,'markerfacecolor','w');
%plotm(43.9,-60,'b+','markersize',10,'linewidth',2)
plotm(36.74,126.13,'bo','markersize',10,'markerfacecolor','b')
cmocean('balance',20);
caxis([-.6 .6])