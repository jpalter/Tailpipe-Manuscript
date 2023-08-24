%carbon_tracker_case_studies_v2.m

%compare fluxes calculated with CT versus other ad hoc methods as published
%in Rutherford et al. for Scotian Shelf (Carioca Buoy, using Sable Island
%Seasonal cycle) and Liu et al. 2022 in the East China Sea.

% Sable Island 43.9N, 59.9 W; %out of tailpipe, source of atm CO2 info
% Carioca buoy 44.3N, 66.3 W;  In Tailpipe

%load  ~/'OneDrive - University of Rhode Island'/Tailpipe/Carbon_tracker_fluxes.mat
%load ~/'OneDrive - University of Rhode Island'/Tailpipe/pCO2_zmxco2_all
%load ~/'Google Drive'/'Shared drives'/'Palter Lab'/Data/'Lanschutzer SOM-FFN'/Patm_PH2O.mat
%load '/Users/jaimepalter/Library/CloudStorage/GoogleDrive-jpalter@uri.edu/Shared drives/Palter Lab/Data/Carbon Tracker/CT2019B/pCO2_zmxco2_all'

indg= ~isnan(pCO2G_all(38,67,:));


%netcdf is in strange format, just import data from text download
%accessed here https://gaw.kishou.go.jp/ on 3/21/2023, saved as
load ~/'OneDrive - University of Rhode Island'/Tailpipe/Sable_island_daily.mat
si_yy = WSA(:,2);
si_mm = WSA(:,3);
si_day = WSA(:,4);
si_xco2 = WSA(:,14);
si_dd = si_yy + si_mm/12 -(1/12) + si_day/365 - (1/365);
%Sable Island 43.9N, 59.9 W
Patm_interp_si = interp1(CTtime(indg),squeeze(Patm_all(41,67,indg)),si_dd);
PH2O_interp_si = interp1(CTtime(indg),squeeze(PH2O_all(41,67,indg)),si_dd);
WAS_pCO2 =  si_xco2 .* (Patm_interp_si - PH2O_interp_si);
figure
subplot(2,1,1);
hold on
plot(si_dd(1:7:end),WAS_pCO2(1:7:end),'b+','linewidth',1)


CTtimeSI = CTtime-1951;
Sable_island_pco2_atm = 273.9904 + 1.9738*CTtimeSI + 5.3132*cos(2*pi*CTtimeSI)...
    + 4.0601*sin(2*pi*CTtimeSI) - 0.5814*cos(4*pi*CTtimeSI) - ...
    2.1740*sin(4*pi*CTtimeSI)-0.5009*cos(6*pi*CTtimeSI)+ 0.5904*sin(6*pi*CTtimeSI);
Sable_island_pco2_atm = Sable_island_pco2_atm';
hold on;
carioca_CTaCO2 = squeeze(pCO2G_all(38,67,:)); %lat/lon 43.5/-67.5
carioca_zm = squeeze(pCO2_zmxco2_all(38,67,:));

indnan = ~isnan(carioca_CTaCO2);
plot(CTtime,carioca_CTaCO2(indnan),'color',[0.5 0.5 0.5])
plot(CTtime,carioca_zm(indnan2),'r')
filt_CTaCO2 = filtfilt(ones(30,1)/30,1,carioca_CTaCO2(indnan));
indnan2 = ~isnan(carioca_zm);
filt_zm = filtfilt(ones(30,1)/30,1,carioca_zm(indnan2));

hold on; plot(CTtime(indnan),filt_CTaCO2,'k','LineWidth',2)
plot(CTtime(indnan2),filt_zm,'r','LineWidth',2)
plot(CTtime,Sable_island_pco2_atm,'b','linewidth',1);
axis([2013 2015 375 430]);
ylabel('pCO_2 (\muAtm)')
grid on

subplot(2,1,2)
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
TAP_xCO2 = table2array(TAP(:,11)); %Mole fraction reported in units of micromol mol-1 (10-6 mol per mol of dry air); abbreviated as ppm (parts per million).
TAP_dd = table2array(TAP(:,2)) + table2array(TAP(:,3))/12 - (1/12) +...
    table2array(TAP(:,4))/365 - (1/365);
Patm_interp = interp1(CTtime,squeeze(Patm_all(102,63,indg)),TAP_dd);
PH2O_interp = interp1(CTtime,squeeze(PH2O_all(102,63,indg)),TAP_dd);
TAP_pCO2 =  TAP_xCO2 .* (Patm_interp - PH2O_interp);


filt_flux_ecs = filtfilt(ones(30,1)/30,1,squeeze(pCO2G_all(102,63,indg)));
filt_flux_ecszm = filtfilt(ones(30,1)/30,1,squeeze(pCO2_zmxco2_all(102,63,indg)));
plot(CTtime',filt_flux_ecs,'k','linewidth',2);
hold on; plot(CTtime,filt_flux_ecszm,'r','linewidth',2);

plot(TAP_dd,TAP_pCO2,'b+','linewidth',1)


axis([2013 2015 375 430]);
ylabel('pCO_2 (\muAtm)')
grid on
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
plotm(36.74,126.13,'bo','markersize',10,'markerfacecolor','w')
cmocean('balance',20);
caxis([-.6 .6])


%calculations for table and response to reviewers
lati(102,63)
loni(102,63)
latil(305,122)
lonil(305,122)
localflux_wCT = squeeze(F_full(305,122,:,:));
localflux_wzm = squeeze(F_xco2zm(305,122,:,:));
this_year = 2014 - 1999;
nanmean(localflux_wCT(:,this_year))
nanmean(localflux_wzm(:,this_year))


%balance of over-estimated uptake in summer when land sink is high and
%underestimate in winter!  Good to address Wanninkhof's comments, but
%annual difference between zonal mean and full is small

%carioca buoy location 38N,67W 
latil(113,134)
lonil(113,134)
localflux_wCT = squeeze(F_full(113,133,:,:));
localflux_wzm = squeeze(F_xco2zm(113,133,:,:));
this_year = 2014 - 1999;
nanmean(localflux_wCT(:,this_year))
nanmean(localflux_wzm(:,this_year))

%calculate RMS error and bias between Sable Island station observations and
%the continuous line drawn fit by Rutherford
SI_interp = interp1(CTtime,Sable_island_pco2_atm,si_dd); %interpolate continuous fit line to time of flask obs
diffe = SI_interp-WAS_pCO2;
RMSd = sqrt(nanmean(diffe.^2))
