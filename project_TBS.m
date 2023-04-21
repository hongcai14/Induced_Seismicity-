%% Hong Ji Cai Term Project

%Note: Numerical values may appear different in comparison to those that
% were presented in class due to minor changes in values due to efficiency
% . Additionally, a few figures will require manual repositioning of their
% respective legend for better visibility. 
%Additionally, ( or ) is non-inclusive while [ or ] is inclusive.
%Also, use Run Section.

%% 1.) Load-In Data/Data Set-Up (Pt. 1)

[~,~,data_1] = ...
    xlsread('Injection_Well_Database_VolumesPressures_OKCOARKNM.xlsx');

[~,~,data_2] = xlsread('The_Human_Induced_Earthquake_Database.xlsx');

%% 2.) Load-In Data/Data Set-Up (Pt. 2) data_1

data_1(1,:) = [];
%^Removes data_1's "Variable Row"

data_1 = cell2mat(data_1);
%^Converts data_1 to numeric array

data_1(data_1(:,3) > -60,:) = [];
%^Constrains data_1 to conterminous USA

data_1(data_1(:,4) < 1932 | data_1(:,4) > 2020 ,:) = [];
%^Constrains data_1 to years (1932,2020); note: min = 1977 and max = 2014

%% 3.) Load-In Data/Data Set-Up (Pt. 3) data_2

data_2 = data_2(strcmp(data_2(:,1),'USA'),:);
%^Constrains data_2 to conterminous USA (no HIEQs (human-induced
%earthquakes) occur in either Arkansas or Hawaii)

%% 4.) Data_1 (Injection Well Database Variables) data_1

longitude_data_1 = data_1(:,3);
latitude_data_1 = data_1(:,2);

injection_year = data_1(:,4);

monthly_injection_rate = mean(data_1(:,6:17),2);
%^units = bbls/month

monthly_injection_pressure = mean(data_1(:,18:29),2);
%^units = psi

data_1(:,6:29) = [];
data_1 = [data_1 monthly_injection_rate monthly_injection_pressure];
%^Replaces all injection rates and injection pressures in data_1 for the 
%injection well's injection year's monthly injection rate and monthly 
%injection pressure 

%% 5.) [Figure 1](Map: Conterminous USA Injection Wells) data_1

figure

m_proj('miller','lon',[-135 -60],'lat',[10 60]);
m_coast('color','k','linewidth',1);
m_grid;

hold on

m_scatter(longitude_data_1,latitude_data_1,'b')

title('Conterminous USA Injections Wells (1932 - 2016)')

%% 6.) data_2 (Induced Earthquake Database Variables) data_2

longitude_data_2 = data_2(:,6);
latitude_data_2 = data_2(:,5);


for ii = 1:length(longitude_data_2)
    if length(longitude_data_2{ii}) == 1
        %^Ignores legitimate geographic coordinate
    else
        longitude_data_2{ii} = longitude_data_2{ii}(2:end);
        latitude_data_2{ii} = latitude_data_2{ii};
        %^Removes '-' in order to use str2num; all latitudes are positive
        %because the US is above the equator
        
        temporary_1 = -str2num(longitude_data_2{ii});
        %^Converts longitude_data_2 to numeric array
        
        longitude_data_2(ii) = {temporary_1};
        %^Re-Converts to cell array
    end
end

index = cell2mat(longitude_data_2) > -60;
%^Identifies geographic locations outside conterminous USA

longitude_data_2(index) = [];
latitude_data_2(index) = [];
data_2(index,:) = [];
%^Constrains longitude_data_2, latitude_data_2, data_2 to conterminous USA

data_2(:,6) = longitude_data_2;
data_2(:,5) = latitude_data_2;
%^Replaces nonuniform geographic coordinates of data_2 with the new,
%uniform ones
    
longitude_data_2 = cell2mat(longitude_data_2);
latitude_data_2 = cell2mat(latitude_data_2);
%^Re-Converts to numeric array 

%% 7.) [Figure 2](Map: Conterminuous HIEQs) data_2

figure

m_proj('miller','lon',[-135 -60],'lat',[10 60]);
m_coast('color','k','linewidth',1);
m_grid;

hold on

m_scatter(longitude_data_2,latitude_data_2)

title('USA Human-Induced Earthquakes (1932 - 2016)')

%% 8.) [Figure 3] Magnitudes vs. HIEQ Causes (Main Class) Scatter Plot

main_cause = data_2(:,2); 
%^HIEQs (human-induced earthquake) causes (main class)

magnitude = cell2mat(data_2(:,13));
%^HIEQs magnitudes

main_cause(isnan(magnitude)) = [];
magnitude(isnan(magnitude)) = [];
%^Removes NaN's corresponding to magnitude

unique(main_cause)
%^Displays unique HIEQ causes (main class)

for ii = 1:length(main_cause)
    if strcmp(main_cause(ii),'CCS')
        main_cause(ii) = {1};
    elseif strcmp(main_cause(ii),'Conventional Oil and Gas')
        main_cause(ii) = {2};
    elseif strcmp(main_cause(ii),'Fracking')
        main_cause(ii) = {3};
    elseif strcmp(main_cause(ii),'Geothermal')
        main_cause(ii) = {4};
    elseif strcmp(main_cause(ii),'Mining')
        main_cause(ii) = {5};
    elseif strcmp(main_cause(ii),'Nuclear explosions')
        main_cause(ii) = {6};
    elseif strcmp(main_cause(ii),'Oil and Gas')
        main_cause(ii) = {7};
    elseif strcmp(main_cause(ii),'Oil and Gas/Waste fluid injection')
        main_cause(ii) = {8};
    elseif strcmp(main_cause(ii),'Research')
        main_cause(ii) = {9};
    elseif strcmp(main_cause(ii),'Waste fluid disposal')
        main_cause(ii) = {10};
    elseif strcmp(main_cause(ii),'Water reservoir impoundment')
        main_cause(ii) = {11};
    end
end
%^Replaces each unique earthquake cause with a respective unique integer
%for plotting

main_cause = cell2mat(main_cause);
%^Converts to numeric array

for ii = 1:11
    if ii < 11
        scatter(main_cause(main_cause == ii),magnitude(...
            main_cause == ii),'filled')
        hold on
    else
        scatter(main_cause(main_cause == ii),magnitude(...
            main_cause == ii),'filled')
    end
end

xlim([0 32])
ylim([-3 8])

xlabel('Earthquake Cause (Main Class)')
ylabel('Magnitude')
title('Magnitude Distribution of HIEQ Causes')

legend('1 = CCS','2 = Conventional Oil and Gas','3 = Fracking',...
    '4 = Geothermal','5 = Mining','6 = Nuclear explosions', ...
    '7 = Oil and Gas','8 = Oil and Gas/Waste fluid injection', ...
    '9 = Research', '10 = Waste fluid disposal',...
    '11 = Water reservoir impoundment')

%% 9.) [Figure 4] Number of Cases vs. Magnitudes (s.t. main causes > 5)

main_cause_frequencies = ones(11,1);
%^Preallocates an array in order to store the frequency of each unique
%earthquake cause

for ii = 1:11
    main_cause_frequencies(ii) = sum(main_cause == ii);
end
%^Calculates the frequency of each unique earthquake cause

copy_main_cause = main_cause;
copy_magnitude = magnitude; 
%^Duplicates the vectors main_cause and magnitude 

select = 1:11;
select = select(main_cause_frequencies < 9);
%^Determines which unique earthquake cause has frequency greater than or
%equal to 9

for ii = select
    copy_magnitude(copy_main_cause == ii) = [];
    copy_main_cause(copy_main_cause == ii) = [];
end
%^Removes values corresponding to frequencies less than 9 from each vector
%Keeps 2,3,5,10,11

select = 1:11;
select = select(main_cause_frequencies > 8);
%^Determines which unique earthquake cause has frequency greater than 8

min(magnitude)
max(magnitude)
%^Determines the minimum and maximum magnitudes

temp_x = -1.5:7.5;
%^Creates vector of x-values

for ii = select
    m = copy_magnitude(copy_main_cause == ii);
    %^Creates vector containing magnitudes corresponding to some unique
    %earthquake cause
    
    temp_y = [sum(m >= -2 & m < -1) sum(m >= -1 & m < 0)...
        sum(m >= 0 & m < 1) sum(m >= 1 & m < 2)...
        sum(m >= 2 & m < 3) sum(m >= 3 & m < 4)...
        sum(m >= 4 & m < 5) sum(m >= 5 & m < 6)...
        sum(m >= 6 & m < 7) sum(m >= 7 & m < 8)];
    %^Calculates frequencies of current unique earthquake cause in
    %intervals [n,n+1), where n ranges from -2 to 7
    
    plot(temp_x,temp_y,'LineWidth',1.5)
    
    if ii < 11
        hold on
    end
end

xlim([-3 13])
ylim([-2 23])

xlabel('Magnitude')
ylabel('Number of Cases')
title('Frequency of Notable HIEQ Causes')

legend('Conventional Oil and Gas','Fracking','Mining',...
    'Waste fluid disposal','Water reservoir impoundment')

%% 10.) [Figure 5] Pie Chart

statistics = ones(11,1);
%^Preallocates an array in order to store the percentage of each unique
%earthquake cause

for ii = 1:11
    A = sum(main_cause == ii)/length(main_cause)*100;
    statistics(ii) = A;
end

pie(statistics)
legend('CCS','Conventional Oil and Gas','Fracking',...
    'Geothermal','Mining','Nuclear explosions', ...
    'Oil and Gas','Oil and Gas/Waste fluid injection', ...
    'Research', 'Waste fluid disposal',...
    'Water reservoir impoundment')

%% 11.) [Figure 6] Temporal Relationship Between Projects and HIEQs

start_year = data_2(:,7);
%^Project (that induced EQ) Start Years

for ii = 1:length(start_year)
    if length(start_year{ii}) == 1
        %^Ignores NaN's and legitimate years
    elseif contains(start_year{ii},',')
        start_year{ii} = NaN;
        %^Replaces double year entries with NaN's, e.g., 1936 , 1956 is
        %replaced with NaN
    elseif strcmp(start_year{ii}(5),'s')
        start_year{ii} = NaN;
        %^Replaces xxxxs with NaN's, e.g., 1900s is replaced with NaN
    elseif strcmp(start_year{ii}(5),' ')
        start_year{ii} = start_year{ii}(1:4);
        %^Removes (month), e.g., 1936 (September) is shortened to 1936
        temporary = str2num(start_year{ii});
        start_year(ii) = {temporary};
    elseif contains(start_year{ii},'/')
        start_year{ii} = start_year{ii}(end-3:end);
        %^Keeps year if format is month/day/year, e.g., 10/24/1996 becomes
        %1996
        temporary = str2num(start_year{ii});
        start_year(ii) = {temporary};
    end
end
%^Corrects irregular data entries   

start_year = cell2mat(start_year);
eq_year = cell2mat(data_2(:,17));

start_year(isnan(eq_year)) = [];
eq_year(isnan(eq_year)) = [];
eq_year(isnan(start_year)) = [];
start_year(isnan(start_year)) = [];
%^Removes NaN's

X = [ones(length(start_year),1) start_year];
%^matrix of Features

[theta,V_theta] = NormalEqn(X,eq_year,[]);
%^Linear Regression (Ordinary Least Squares)

scatter(start_year,eq_year,'b')
%^Plots Data Points

hold on

p1 = plot([1900 2020],theta(1) + theta(2)*[1900 2020],'r','LineWidth',...
    1.25);
%Plots Regression Function

hold on

p2 = plot([1900 2020],[1900 2020],'LineStyle','--','Color','[1 0.5 0]'...
    ,'LineWidth',1.25);
%^Plots Identity Function; if a Data Point lies below this line, then
%something is amiss

ylim([1900 2020])

xlabel('Project Start Year')
ylabel('HIEQ Year')
title('Temporal Relationship Between Projects and HIEQs')

legend('Data Points', 'Regression Function', 'Identity Function')

%% 12.) Injection Rates and Number of Induced Earthquakes

slide_3_pic_1 = cell2mat(data_2(:,[5 6 13 17]));
%^Creates an abridged array of data_2: longitude, latitude, magnitude, year

slide_3_pic_1(slide_3_pic_1(:,2) < min(longitude_data_1) | ...
    slide_3_pic_1(:,2) > max(longitude_data_1),:) = [];
slide_3_pic_1(slide_3_pic_1(:,1) < min(latitude_data_1) | ...
    slide_3_pic_1(:,1) > max(latitude_data_1),:) = [];
%^Constrains HIEQs to quadrangle enclosed by data_1's minimum longitude,
%minimum latitude, maximum longitude, maximum latitude

%% 13.) [Figure 7] Slide 3 Part 1 (Injection Rate)(Raw) 

slide_3_pic_1(isnan(slide_3_pic_1(:,4)),:) = [];
%^Removes rows containing NaN(s)

data_1(data_1(:,6) <= 0,:) = [];
%^Removes rows from data_1 containing injection rates less than or equal to
%0

data_1 = sortrows(data_1,4);
%^Sorts rows in data_1 with respect to years (column 4) from 1977 to 2014

r1975_1980 = [0 0];
r1980_1985 = [0 0];
r1985_1990 = [0 0];
r1990_1995 = [0 0];
r1995_2000 = [0 0];
r2000_2005 = [0 0];
r2005_2010 = [0 0];
r2010_2015 = [0 0];
%^Preallocates space, rx_y = [x,y), where 
%   x = the sum of injection rates
%   y = the number of injection wells

for ii = 1:length(data_1)
    if data_1(ii,4) >= 1975 && data_1(ii,4) < 1980
        r1975_1980 = r1975_1980 + [data_1(ii,6) 1];
    elseif data_1(ii,4) >= 1980 && data_1(ii,4) < 1985
        r1980_1985 = r1980_1985 + [data_1(ii,6) 1];
    elseif data_1(ii,4) >= 1985 && data_1(ii,4) < 1990
        r1985_1990 = r1985_1990 + [data_1(ii,6) 1];
    elseif data_1(ii,4) >= 1990 && data_1(ii,4) < 1995
        r1990_1995 = r1990_1995 + [data_1(ii,6) 1];
    elseif data_1(ii,4) >= 1995 && data_1(ii,4) < 2000
        r1995_2000 = r1995_2000 + [data_1(ii,6) 1];
    elseif data_1(ii,4) >= 2000 && data_1(ii,4) < 2005
        r2000_2005 = r2000_2005 + [data_1(ii,6) 1];
    elseif data_1(ii,4) >= 2005 && data_1(ii,4) < 2010
        r2005_2010 = r2005_2010 + [data_1(ii,6) 1];
    else
        r2010_2015 = r2010_2015 + [data_1(ii,6) 1];
    end
end
%^Calculates the sum of injection rates and the number of injection wells
%from [1975 + n,1980 + n), where n = 0:5:35

injection_rates_mean = [r1975_1980(1)/r1975_1980(2) ...
    r1980_1985(1)/r1980_1985(2) r1985_1990(1)/r1985_1990(2) ...
    r1990_1995(1)/r1990_1995(2) r1995_2000(1)/r1995_2000(2) ...
    r2000_2005(1)/r2000_2005(2) r2005_2010(1)/r2005_2010(2) ...
    r2010_2015(1)/r2010_2015(2)];
%^Calculates the average injection rate for each interval I_n = [1975 +
%n,1980 + n), where n = 0:5:35.

slide_3_pic_1 = sortrows(slide_3_pic_1,4);
%^Sorts rows in slide_3_pic_1 with respect to years (column 4) from 1952 to
%2016

d_slide_3_pic_1 = slide_3_pic_1;
%^Duplicates slide_3_pic_1

EQ_1975_1980 = [-10 0];
EQ_1980_1985 = [-10 0];
EQ_1985_1990 = [-10 0];
EQ_1990_1995 = [-10 0];
EQ_1995_2000 = [-10 0];
EQ_2000_2005 = [-10 0];
EQ_2005_2010 = [-10 0];
EQ_2010_2015 = [-10 0];
%^Preallocates space EQ_x_y, where
%   x = maximum magnitude from the five-year half-closed interval [x,y)
%   y = the number of earthquakes from the five-year half-closed interval
%   [x,y)

for ii = 1:length(d_slide_3_pic_1)
    if d_slide_3_pic_1(ii,4) >= 1975 && d_slide_3_pic_1(ii,4) < 1980
        EQ_1975_1980 = EQ_1975_1980 + [0 1];
        if d_slide_3_pic_1(ii,3) > EQ_1975_1980(1)
            EQ_1975_1980(1) = d_slide_3_pic_1(ii,3);
        end
    elseif d_slide_3_pic_1(ii,4) >= 1980 && d_slide_3_pic_1(ii,4) < 1985
        EQ_1980_1985 = EQ_1980_1985 + [0 1];
        if d_slide_3_pic_1(ii,3) > EQ_1980_1985(1)
            EQ_1980_1985(1) = d_slide_3_pic_1(ii,3);
        end
    elseif d_slide_3_pic_1(ii,4) >= 1985 && d_slide_3_pic_1(ii,4) < 1990
        EQ_1985_1990 = EQ_1985_1990 + [0 1];
        if d_slide_3_pic_1(ii,3) > EQ_1985_1990(1)
            EQ_1985_1990(1) = d_slide_3_pic_1(ii,3);
        end
    elseif d_slide_3_pic_1(ii,4) >= 1990 && d_slide_3_pic_1(ii,4) < 1995
        EQ_1990_1995 = EQ_1990_1995 + [0 1];
        if d_slide_3_pic_1(ii,3) > EQ_1990_1995(1)
            EQ_1990_1995(1) = d_slide_3_pic_1(ii,3);
        end
    elseif d_slide_3_pic_1(ii,4) >= 1995 && d_slide_3_pic_1(ii,4) < 2000
        EQ_1995_2000 = EQ_1995_2000 + [0 1];
        if d_slide_3_pic_1(ii,3) > EQ_1995_2000(1)
            EQ_1995_2000(1) = d_slide_3_pic_1(ii,3);
        end
    elseif d_slide_3_pic_1(ii,4) >= 2000 && d_slide_3_pic_1(ii,4) < 2005
        EQ_2000_2005 = EQ_2000_2005 + [0 1];
        if d_slide_3_pic_1(ii,3) > EQ_2000_2005(1)
            EQ_2000_2005(1) = d_slide_3_pic_1(ii,3);
        end
    elseif d_slide_3_pic_1(ii,4) >= 2005 && d_slide_3_pic_1(ii,4) < 2010
        EQ_2005_2010 = EQ_2005_2010 + [0 1];
        if d_slide_3_pic_1(ii,3) > EQ_2005_2010(1)
            EQ_2005_2010(1) = d_slide_3_pic_1(ii,3);
        end
    elseif d_slide_3_pic_1(ii,4) >= 2010 && d_slide_3_pic_1(ii,4) < 2015    
        EQ_2010_2015 = EQ_2010_2015 + [0 1];
        if d_slide_3_pic_1(ii,3) > EQ_2010_2015(1)
            EQ_2010_2015(1) = d_slide_3_pic_1(ii,3);
        end
    end
end
%^Determines the maximum magnitude and the number of earthquakes for each
%five-year half-closed interval 

EQ_1975_2015 = [EQ_1975_1980' EQ_1980_1985' EQ_1985_1990'...
                EQ_1990_1995' EQ_1995_2000' EQ_2000_2005'...
                EQ_2005_2010' EQ_2010_2015'];
            
injection_rates_mean(:,EQ_1975_2015(2,:) == 0) = [];
EQ_1975_2015(:,EQ_1975_2015(2,:) == 0) = [];
%^Removes columns corresponding to five-year half-closed intervals
%unassociated with any HIEQs
            
EQ_1975_2015(:,isnan(injection_rates_mean)) = [];
injection_rates_mean(:,isnan(injection_rates_mean)) = [];
%^Removes columns containing NaN's 

maximum_magnitudes_1975_2015 = EQ_1975_2015(1,:);
number_of_earthquakes_1975_2015 = EQ_1975_2015(2,:);

yyaxis left
plot(injection_rates_mean,maximum_magnitudes_1975_2015)
ylabel('Maximum "Observed" Magnitude')
xlabel('Mean Injection Rate (bbls per month)')
ylim([0 14])

yyaxis right
plot(injection_rates_mean,number_of_earthquakes_1975_2015)
ylabel('Number of Earthquakes')

title('Injection Rate Effects')

%% 14.) [Figure 8] Slide 3 Part 2 (Injection Rate)(Linear Regression)

irm = injection_rates_mean';
%^Duplicates injection_rates_mean

X_1 = [ones(5,1) irm irm.^2];
%^matrix of Features

[theta_1,V_theta_1] = NormalEqn(X_1,number_of_earthquakes_1975_2015',[]);
%^Linear Regression (Ordinary Least Squares), the numbers of earthquakes

syms c

c = solve(c + theta_1(2)*irm(3) + theta_1(3)*irm(3)^2 == 1,c);
%^Solves for c, the y-intercept, since theta_1(0) fails due to the
%piecewise approach 

X_2 = [X_1 irm.^3 irm.^4];
%^matrix of Features

[theta_2,V_theta_2] = NormalEqn(X_2,maximum_magnitudes_1975_2015',[]);
%^Linear Regression (Ordinary Least Squares), the maximum magnitudes

x_1 = irm(3):100:2.4*10^4;
%^Creates a vector of x values for number_of_earthquakes

x_2 = 0.8*10^4:100:2.4*10^4;
%^Creates a vector of x values for maximum_magnitudes

yyaxis right

plot(x_1,c + theta_1(2)*x_1 + theta_1(3)*x_1.^2,...
    'LineStyle','-','Color','[1 0.5 0]')

hold on

plot([0 irm(3)],[1 1],...
    'LineStyle','-','Color','[1 0.5 0]')

xlim([0.8*10^4 2.4*10^4])
ylim([0 14])

xlabel('Mean Injection Rate (bbls per month)')
ylabel('Number of Earthquakes')
title('Injection Rate Effects')

legend('f(x) = 1.4481e-07x^2 - 0.0039x + 24.1600')

yyaxis left

plot(x_2, theta_2(1) + theta_2(2)*x_2 + theta_2(3)*x_2.^2 ...
    + theta_2(4)*x_2.^3 + theta_2(5)*x_2.^4,...
    'LineStyle','-','Color','b')

ylim([0 14])

ylabel('Maximum "Observed" Magnitude')

%% 15.) [Figure 9] Slide 3 Part 3 (Injection Pressure)(Raw)

%Note: I ignore the five-year half-closed intervals [1980,1985),
%[1990,1995), [1995,2000) since 

p1975_1980 = [0 0];
p1985_1990 = [0 0];
p2000_2005 = [0 0];
p2005_2010 = [0 0];
p2010_2015 = [0 0];
%^Preallocates space, px_y = [x,y), where 
%   x = the sum of injection pressures
%   y = the number of injection wells

for ii = 1:length(data_1)
    if data_1(ii,4) >= 1975 && data_1(ii,4) < 1980
        p1975_1980 = p1975_1980 + [data_1(ii,7) 1];
    elseif data_1(ii,4) >= 1985 && data_1(ii,4) < 1990
        p1985_1990 = p1985_1990 + [data_1(ii,7) 1];
    elseif data_1(ii,4) >= 2000 && data_1(ii,4) < 2005
        p2000_2005 = p2000_2005 + [data_1(ii,7) 1];
    elseif data_1(ii,4) >= 2005 && data_1(ii,4) < 2010
        p2005_2010 = p2005_2010 + [data_1(ii,7) 1];
    elseif data_1(ii,4) >= 2010 && data_1(ii,4) < 2015
        p2010_2015 = p2010_2015 + [data_1(ii,7) 1];
    end
end
%^Calculates the sum of injection pressures and the number of injection 
%wells from [1975 + n,1980 + n), where n = 0:5:35

injection_pressures_mean = [p1975_1980(1)/p1975_1980(2) ...
    p1985_1990(1)/p1985_1990(2) p2000_2005(1)/p2000_2005(2)...
    p2005_2010(1)/p2005_2010(2) p2010_2015(1)/p2010_2015(2)]';
%^Calculates the average injection pressure for each interval I_n = [1975 +
%n,1980 + n), where n = 0:5:35.

pressure_matrix = [injection_pressures_mean...
    maximum_magnitudes_1975_2015' number_of_earthquakes_1975_2015'];
%^Stores relevant information for subsequent plot

pressure_matrix = sortrows(pressure_matrix,1);
%^Sorts rows in pressure_matrix with respect to average injection pressures
%(column 1)

yyaxis left

plot(pressure_matrix(:,1),pressure_matrix(:,2))

ylim([0 14])

ylabel('Maximum "Observed" Magnitude')
xlabel('Mean Injection Pressure (psi)')

yyaxis right

plot(pressure_matrix(:,1),pressure_matrix(:,3))

ylabel('Number of Earthquakes')
title('Injection Pressure Effects')

%% 16.) Slide 3 Part 4 (Injection Pressure)(Linear Regression)

%Note: I included the point (600,10) to my linear regression of the number
%of earthquakes vs. injection pressure plot, which I manually selected from
%the previous, discrete plot in order to form a more accurate "fit"

ipm = [injection_pressures_mean; 600];
%^Duplicates injection_pressures_mean; adds x-coordinate of additional
%point (600,10)

noe_1975_2015 = [number_of_earthquakes_1975_2015 10]';
%^Duplicates number_of_earthquakes_1975_2015; adds y-coordinate of
%additional point (600,10)

X_3 = [ones(4,1) ipm(3:6) ipm(3:6).^2];
%^matrix of Features
    
[theta_3,V_theta_3] = NormalEqn(X_3,noe_1975_2015(3:6),[]);
%^Linear Regression (Ordinary Least Squares), the numbers of earthquakes

syms k

k = solve(k + theta_3(2)*ipm(3) + theta_3(3)*ipm(3)^2 == 1,k);
%^Solves for k, the y-intercept, since theta_3(0) fails due to the
%piecewise approach 

X_4 = [ones(5,1) ipm(1:5) ipm(1:5).^2 ipm(1:5).^3];
%^matrix of Features

[theta_4,V_theta_4] = NormalEqn(X_4,maximum_magnitudes_1975_2015',[]);
%^Linear Regression (Ordinary Least Squares), the maximum magnitudes

x_3 = ipm(3):1:ipm(5);
%^Creates a vector of x values for number_of_earthquakes

x_4 = 0:700;
%^Creates a vector of x values for maximum_magnitudes

yyaxis right

plot(x_3,k + theta_3(2)*x_3 + theta_3(3)*x_3.^2,...
    'LineStyle','-','Color','[1 0.5 0]')

hold on

plot([0 ipm(3)],[1 1],'LineStyle','-','Color','[1 0.5 0]')

legend('f(x) = -5.0694e-05x^2 + 0.1077x - 36.0613')

xlim([0 700])
ylim([0 14])

xlabel('Mean Injection Pressure (psi)')
ylabel('Number of Earthquakes')
title('Injection Pressure Effects')

yyaxis left

plot(x_4, theta_4(1) + theta_4(2)*x_4 + theta_4(3)*x_4.^2 ...
+ theta_4(4)*x_4.^3,'LineStyle','-','Color','b')

ylim([0 14])

ylabel('Maximum "Observed" Magnitude')
