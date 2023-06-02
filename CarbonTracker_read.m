% Read in carbon tracker (CT2019B)
% 1x1 degree, time resolution 3 hours
% ECMWF operational model for wind
% gas transfer velocity from Wanninkhof 1992
% Transport Model 5 for atm transport of CO2
load /Users/sarahnickford/Documents/Saildrone/processed/Saildrone2019_pH_pCO2_Final lon_10min_ph lat_10min_ph mtime_10minpH pCO2_atm_10min lon_asvco2 lat_asvco2 mtime_asvco2_trim

cd /Users/sarahnickford/Documents/Saildrone/CO2/CT2019B

% contents:
% air_mass = air mass (kg)
% blh = planetary boundary layer thickness (m)
% boundary = boundary (no unit, positive up)
% calendar_components = (no unit)
% co2 = mole fraction of carbon dioxide in air (umol/mol)- CT2019B estimate of total atmospheric carbon dioxide
% decimal_date = decimal date (years)
% gph = geopotential height at level boundaries (m)
% latitude (degrees north)
% level (no unit, positive up), level 1 ~ 25m
% longitude (degrees east)
% orography = surface geopotential (m^2/s^2)
% pressure = air pressure at level boundaries (Pa)
% specific_humidity = mass fraction of water vapor in moist air; (kg water)/(kg air)
% temperature = air temperature at level center (K)
% time = (days since 2000-01-01 00:00:00 UTC)
% time_components = year, month, day, hour, minute, second
% u = zonal wind (m/s)
% v = meridional wind (m/s)

%% Read in CT2019B North Atlantic ocean files
files = dir('*.nc');
for ii = 1:length(files)
    ct2019b(ii) = netcdf_load(files(ii).name);
end

% convert time to mtime
for ik = 1:length(files)
    ct2019b(ik).mtime = double(datenum('2000', 'yyyy') + ct2019b(ik).time);% / 864e2;
end

for io = 1:length(files)
    if io == 1
%         longitude = ct2019b(io).longitude;
%         latitude = ct2019b(io).latitude;
        mtime = ct2019b(io).mtime;
        surf_xco2 = double(permute(squeeze(ct2019b(io).co2(:,:,1,:)),[2 1 3]));
        pressure = double(permute(squeeze(ct2019b(io).pressure(:,:,1,:)),[2 1 3]));
        specific_humidity = double(permute(squeeze(ct2019b(io).specific_humidity(:,:,1,:)),[2 1 3]));
    else
        this_surf_xco2 = double(permute(squeeze(ct2019b(io).co2(:,:,1,:)),[2 1 3]));
        this_surf_press = double(permute(squeeze(ct2019b(io).pressure(:,:,1,:)),[2 1 3]));
        this_surf_humid = double(permute(squeeze(ct2019b(io).specific_humidity(:,:,1,:)),[2 1 3]));
%         longitude = cat(1,longitude,ct2019b(io).longitude);
%         latitude = cat(1,latitude,ct2019b(io).latitude);
        mtime = cat(1,mtime,ct2019b(io).mtime);
        surf_xco2 = cat(3,surf_xco2,this_surf_xco2);
        pressure = cat(3,pressure,this_surf_press);
        specific_humidity = cat(3,specific_humidity,this_surf_humid);
    end
end
longitude = ct2019b(1).longitude;
latitude = ct2019b(1).latitude;

surf_xco2_interp = interp3(longitude,latitude,mtime,surf_xco2,lon_asvco2,lat_asvco2,mtime_asvco2_trim);
% need to convert to pCO2
% adjust from 25m to 0m? Possibly unnecessary as Tudor Hill tower is higher

%% Read in Tudor Hill observations
file1='/Users/sarahnickford/Documents/Saildrone/CO2/TudorHill_CO2/co2_bmw_surface-flask_1_ccgg_event2020.xlsx';
idxtimeval = [datenum(2019,1,20,0,0,0) datenum(2019,3,05,0,0,0)];
tudorhill = importdata(file1);
tudorhill.atm_xco2 = tudorhill.data(:,11);
tudorhill.atm_xco2_uncertainty = tudorhill.data(:,12);
tudorhill.mtime = datenum(tudorhill.data(:,1),tudorhill.data(:,2),...
        tudorhill.data(:,3),tudorhill.data(:,4),tudorhill.data(:,5),00);
atm_xco2_SD_idx = find(tudorhill.mtime >= idxtimeval(1) & tudorhill.mtime <= idxtimeval(2));
atm_xco2_TH_SD = tudorhill.atm_xco2(atm_xco2_SD_idx);
atm_xco2_TH_SDuncert = tudorhill.atm_xco2_uncertainty(atm_xco2_SD_idx);
co2_time = tudorhill.mtime(atm_xco2_SD_idx);
co2_time = co2_time(1:2:end);
atm_xco2_TH_SD = atm_xco2_TH_SD(1:2:end);
atm_xco2_TH_SDuncert = atm_xco2_TH_SDuncert(1:2:end);

%% read in the 3x2 global file
cd /Users/sarahnickford/Documents/Saildrone/CO2/CT-NRT/CTNRT_2019/
files = dir('*.nc');
for ii = 1:length(files)
    ctnrt(ii) = netcdf_load(files(ii).name);
end

% convert time to mtime
for ik = 1:length(files)
    ctnrt(ik).mtime = double(datenum('2000', 'yyyy') + ctnrt(ik).time);% / 864e2;
end

for io = 1:length(files)
    if io == 1
%         longitude = ct2019b(io).longitude;
%         latitude = ct2019b(io).latitude;
        ctnrt_mtime = ctnrt(io).mtime;
        ctnrt_surf_xco2 = double(permute(squeeze(ctnrt(io).co2(:,:,1,:)),[2 1 3]));
        ctnrt_pressure = double(permute(squeeze(ctnrt(io).pressure(:,:,1,:)),[2 1 3]));
        ctnrt_specific_humidity = double(permute(squeeze(ctnrt(io).specific_humidity(:,:,1,:)),[2 1 3]));
    else
        this_surf_xco2 = double(permute(squeeze(ctnrt(io).co2(:,:,1,:)),[2 1 3]));
        this_surf_press = double(permute(squeeze(ctnrt(io).pressure(:,:,1,:)),[2 1 3]));
        this_surf_humid = double(permute(squeeze(ctnrt(io).specific_humidity(:,:,1,:)),[2 1 3]));
%         longitude = cat(1,longitude,ct2019b(io).longitude);
%         latitude = cat(1,latitude,ct2019b(io).latitude);
        ctnrt_mtime = cat(1,ctnrt_mtime,ctnrt(io).mtime);
        ctnrt_surf_xco2 = cat(3,ctnrt_surf_xco2,this_surf_xco2);
        ctnrt_pressure = cat(3,ctnrt_pressure,this_surf_press);
        ctnrt_specific_humidity = cat(3,ctnrt_specific_humidity,this_surf_humid);
    end
end
ctnrt_longitude = ctnrt(1).longitude;
ctnrt_latitude = ctnrt(1).latitude;

ctnrt_surf_xco2_interp = interp3(ctnrt_longitude,ctnrt_latitude,ctnrt_mtime,ctnrt_surf_xco2,lon_asvco2,lat_asvco2,mtime_asvco2_trim);

%% Plot
figure;hold on
plot(mtime_asvco2_trim,xco2_air_asvco2_trim,'k','linewidth',2)
plot(mtime_asvco2_trim,surf_xco2_interp,'linewidth',2)
plot(mtime_asvco2_trim,ctnrt_surf_xco2_interp,'linewidth',2)
plot(co2_time,atm_xco2_TH_SD,'.m','markersize',45)
datetick('x','mm/dd','keepticks')
legend('Saildrone xCO2','CT2019b 1x1nam xCO2','CT-NRT 3x2glo','Tudor Hill xCO2')
set(gca,'fontsize',25);
ylabel('xCO2 [\mumol/mol]')
RMSerr = sqrt((sum((surf_xco2_interp-xco2_air_asvco2_trim).^2,'omitnan'))/size(xco2_air_asvco2_trim,1));
% 1.9 umol/mol
RMSerr = sqrt((sum((ctnrt_surf_xco2_interp-xco2_air_asvco2_trim).^2,'omitnan'))/size(xco2_air_asvco2_trim,1));
% 2.2 umol/mol



% save /Users/sarahnickford/Documents/Saildrone/processed/SD2019_atmCO2_products ...
%     ctnrt ctnrt_surf_xco2_interp ctnrt_longitude ctnrt_latitude co2_time ...
%     atm_xco2_TH_SD atm_xco2_TH_SDuncert surf_xco2_interp latitude longitude ...
%     mtime specific_humidity surf_xco2 pressure ct2019b tudorhill -v7.3

%% Level
% model level || Mean height (m)
% 1	                25	
% 2	                103	
% 3	                247	
% 4	                480	
% 5	                814	
% 6	                1259	
% 7	                1822	
% 8	                2508	
% 9	                3317	
% 10	            4248	
% 11	            5300	
% 12	            6467	
% 13	            7741	
% 14	            9114
% 15	            10588
% 16	            12184
% 17	            13928
% 18	            15843
% 19	            17983
% 20	            20412
% 21	            24433
% 22	            30003
% 23	            35895
% 24	            43210
% 25	            123622