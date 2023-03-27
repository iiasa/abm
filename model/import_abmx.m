%% Import data of ABMX

clear;
year_=2010;
number_years=10;
number_quarters=4*number_years;
quarters_num=datenum(year_,4:3:(number_years+1)*12+1,0);
horizon=12;
number_seeds=500;
number_sectors=62;

% Load data
load('../data/data/1996Q1-2019Q4.mat');

for i=1:number_quarters
    quarter_num=quarters_num(i);
    year_num=datenum(datetime(year(datetime(datestr(quarter_num)))+1,1,0));
    
    % Load data
    load(['1_',num2str(year(datetime(datestr(quarters_num(i))))),'Q',num2str(quarter(datetime(datestr(quarters_num(i))))),'.mat']);
    
    forecast_quarter_num=datenum(datetime(year(datetime(datestr(quarter_num))),month(datetime(datestr(quarter_num)))+(3*horizon)+1,0));
    q=quarter(datetime(datestr(quarters_num(i))));
    
    real_gdp_growth_quarterly=diff(log(real_gdp));
    real_gdp_growth_quarterly=exp(real_gdp_growth_quarterly)-1;
    real_gdp_quarterly=data.real_gdp_quarterly(data.quarters_num==quarter_num)*cumprod(1+real_gdp_growth_quarterly);
    
    model.real_gdp_quarterly=[repmat(data.real_gdp_quarterly(data.quarters_num==quarter_num),1,number_seeds);real_gdp_quarterly];
    model.real_gdp=[repmat(data.real_gdp(data.years_num==year_num),1,number_seeds);toannual(real_gdp_quarterly(5-q:end-mod(q,4),:)')'];
    model.real_gdp_growth=diff(log([repmat(data.real_gdp(data.years_num==year_num),1,number_seeds);toannual(real_gdp_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.real_gdp_growth=exp(model.real_gdp_growth)-1;
    model.real_gdp_growth=[repmat(data.real_gdp_growth(data.years_num==year_num),1,number_seeds);model.real_gdp_growth];
    
    model.real_gdp_growth_quarterly=diff(log([repmat(data.real_gdp_quarterly(data.quarters_num==quarter_num),1,number_seeds);real_gdp_quarterly]));
    % calculate discrete compounding rate
    model.real_gdp_growth_quarterly=exp(model.real_gdp_growth_quarterly)-1;
    model.real_gdp_growth_quarterly=[repmat(data.real_gdp_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.real_gdp_growth_quarterly];
    
    nominal_gdp_growth_quarterly=diff(log(nominal_gdp));
    nominal_gdp_growth_quarterly=exp(nominal_gdp_growth_quarterly)-1;
    nominal_gdp_quarterly=data.nominal_gdp_quarterly(data.quarters_num==quarter_num)*cumprod(1+nominal_gdp_growth_quarterly);
    
    model.nominal_gdp_quarterly=[repmat(data.nominal_gdp_quarterly(data.quarters_num==quarter_num),1,number_seeds);nominal_gdp_quarterly];
    model.nominal_gdp=[repmat(data.nominal_gdp(data.years_num==year_num),1,number_seeds);toannual(nominal_gdp_quarterly(5-q:end-mod(q,4),:)')'];
    model.nominal_gdp_growth=diff(log([repmat(data.nominal_gdp(data.years_num==year_num),1,number_seeds);toannual(nominal_gdp_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.nominal_gdp_growth=exp(model.nominal_gdp_growth)-1;
    model.nominal_gdp_growth=[repmat(data.nominal_gdp_growth(data.years_num==year_num),1,number_seeds);model.nominal_gdp_growth];
    
    model.nominal_gdp_growth_quarterly=diff(log([repmat(data.nominal_gdp_quarterly(data.quarters_num==quarter_num),1,number_seeds);nominal_gdp_quarterly]));
    % calculate discrete compounding rate
    model.nominal_gdp_growth_quarterly=exp(model.nominal_gdp_growth_quarterly)-1;
    model.nominal_gdp_growth_quarterly=[repmat(data.nominal_gdp_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.nominal_gdp_growth_quarterly];
    
    gdp_deflator_quarterly=nominal_gdp_quarterly./real_gdp_quarterly;
      
    model.gdp_deflator_quarterly=[repmat(data.gdp_deflator_quarterly(data.quarters_num==quarter_num),1,number_seeds);gdp_deflator_quarterly];
    model.gdp_deflator=[repmat(data.gdp_deflator(data.years_num==year_num),1,number_seeds);toannual_mean(gdp_deflator_quarterly(5-q:end-mod(q,4),:)')'];
    model.gdp_deflator_growth=diff(log([repmat(data.gdp_deflator(data.years_num==year_num),1,number_seeds);toannual_mean(gdp_deflator_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.gdp_deflator_growth=exp(model.gdp_deflator_growth)-1;
    model.gdp_deflator_growth=[repmat(data.gdp_deflator_growth(data.years_num==year_num),1,number_seeds);model.gdp_deflator_growth];
    
    model.gdp_deflator_growth_quarterly=diff(log([repmat(data.gdp_deflator_quarterly(data.quarters_num==quarter_num),1,number_seeds);gdp_deflator_quarterly]));
    % calculate discrete compounding rate
    model.gdp_deflator_growth_quarterly=exp(model.gdp_deflator_growth_quarterly)-1;
    model.gdp_deflator_growth_quarterly=[repmat(data.gdp_deflator_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.gdp_deflator_growth_quarterly];
    
    real_gva_growth_quarterly=diff(log(real_gva));
    real_gva_growth_quarterly=exp(real_gva_growth_quarterly)-1;
    real_gva_quarterly=data.real_gva_quarterly(data.quarters_num==quarter_num)*cumprod(1+real_gva_growth_quarterly);
    
    model.real_gva_quarterly=[repmat(data.real_gva_quarterly(data.quarters_num==quarter_num),1,number_seeds);real_gva_quarterly];
    model.real_gva=[repmat(data.real_gva(data.years_num==year_num),1,number_seeds);toannual(real_gva_quarterly(5-q:end-mod(q,4),:)')'];
    model.real_gva_growth=diff(log([repmat(data.real_gva(data.years_num==year_num),1,number_seeds);toannual(real_gva_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.real_gva_growth=exp(model.real_gva_growth)-1;
    model.real_gva_growth=[repmat(data.real_gva_growth(data.years_num==year_num),1,number_seeds);model.real_gva_growth];
    
    model.real_gva_growth_quarterly=diff(log([repmat(data.real_gva_quarterly(data.quarters_num==quarter_num),1,number_seeds);real_gva_quarterly]));
    % calculate discrete compounding rate
    model.real_gva_growth_quarterly=exp(model.real_gva_growth_quarterly)-1;
    model.real_gva_growth_quarterly=[repmat(data.real_gva_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.real_gva_growth_quarterly];
    
    nominal_gva_growth_quarterly=diff(log(nominal_gva));
    nominal_gva_growth_quarterly=exp(nominal_gva_growth_quarterly)-1;
    nominal_gva_quarterly=data.nominal_gva_quarterly(data.quarters_num==quarter_num)*cumprod(1+nominal_gva_growth_quarterly);
    
    model.nominal_gva_quarterly=[repmat(data.nominal_gva_quarterly(data.quarters_num==quarter_num),1,number_seeds);nominal_gva_quarterly];
    model.nominal_gva=[repmat(data.nominal_gva(data.years_num==year_num),1,number_seeds);toannual(nominal_gva_quarterly(5-q:end-mod(q,4),:)')'];
    model.nominal_gva_growth=diff(log([repmat(data.nominal_gva(data.years_num==year_num),1,number_seeds);toannual(nominal_gva_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.nominal_gva_growth=exp(model.nominal_gva_growth)-1;
    model.nominal_gva_growth=[repmat(data.nominal_gva_growth(data.years_num==year_num),1,number_seeds);model.nominal_gva_growth];
    
    model.nominal_gva_growth_quarterly=diff(log([repmat(data.nominal_gva_quarterly(data.quarters_num==quarter_num),1,number_seeds);nominal_gva_quarterly]));
    % calculate discrete compounding rate
    model.nominal_gva_growth_quarterly=exp(model.nominal_gva_growth_quarterly)-1;
    model.nominal_gva_growth_quarterly=[repmat(data.nominal_gva_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.nominal_gva_growth_quarterly];
    
    gva_deflator_quarterly=nominal_gva_quarterly./real_gva_quarterly;
      
    model.gva_deflator_quarterly=[repmat(data.gva_deflator_quarterly(data.quarters_num==quarter_num),1,number_seeds);gva_deflator_quarterly];
    model.gva_deflator=[repmat(data.gva_deflator(data.years_num==year_num),1,number_seeds);toannual_mean(gva_deflator_quarterly(5-q:end-mod(q,4),:)')'];
    model.gva_deflator_growth=diff(log([repmat(data.gva_deflator(data.years_num==year_num),1,number_seeds);toannual_mean(gva_deflator_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.gva_deflator_growth=exp(model.gva_deflator_growth)-1;
    model.gva_deflator_growth=[repmat(data.gva_deflator_growth(data.years_num==year_num),1,number_seeds);model.gva_deflator_growth];
    
    model.gva_deflator_growth_quarterly=diff(log([repmat(data.gva_deflator_quarterly(data.quarters_num==quarter_num),1,number_seeds);gva_deflator_quarterly]));
    % calculate discrete compounding rate
    model.gva_deflator_growth_quarterly=exp(model.gva_deflator_growth_quarterly)-1;
    model.gva_deflator_growth_quarterly=[repmat(data.gva_deflator_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.gva_deflator_growth_quarterly];
    
    real_household_consumption_growth_quarterly=diff(log(real_household_consumption));
    real_household_consumption_growth_quarterly=exp(real_household_consumption_growth_quarterly)-1;
    real_household_consumption_quarterly=data.real_household_consumption_quarterly(data.quarters_num==quarter_num)*cumprod(1+real_household_consumption_growth_quarterly);
    
    model.real_household_consumption_quarterly=[repmat(data.real_household_consumption_quarterly(data.quarters_num==quarter_num),1,number_seeds);real_household_consumption_quarterly];
    model.real_household_consumption=[repmat(data.real_household_consumption(data.years_num==year_num),1,number_seeds);toannual(real_household_consumption_quarterly(5-q:end-mod(q,4),:)')'];
    model.real_household_consumption_growth=diff(log([repmat(data.real_household_consumption(data.years_num==year_num),1,number_seeds);toannual(real_household_consumption_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.real_household_consumption_growth=exp(model.real_household_consumption_growth)-1;
    model.real_household_consumption_growth=[repmat(data.real_household_consumption_growth(data.years_num==year_num),1,number_seeds);model.real_household_consumption_growth];
    
    model.real_household_consumption_growth_quarterly=diff(log([repmat(data.real_household_consumption_quarterly(data.quarters_num==quarter_num),1,number_seeds);real_household_consumption_quarterly]));
    % calculate discrete compounding rate
    model.real_household_consumption_growth_quarterly=exp(model.real_household_consumption_growth_quarterly)-1;
    model.real_household_consumption_growth_quarterly=[repmat(data.real_household_consumption_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.real_household_consumption_growth_quarterly];
    
    nominal_household_consumption_growth_quarterly=diff(log(nominal_household_consumption));
    nominal_household_consumption_growth_quarterly=exp(nominal_household_consumption_growth_quarterly)-1;
    nominal_household_consumption_quarterly=data.nominal_household_consumption_quarterly(data.quarters_num==quarter_num)*cumprod(1+nominal_household_consumption_growth_quarterly);
    
    model.nominal_household_consumption_quarterly=[repmat(data.nominal_household_consumption_quarterly(data.quarters_num==quarter_num),1,number_seeds);nominal_household_consumption_quarterly];
    model.nominal_household_consumption=[repmat(data.nominal_household_consumption(data.years_num==year_num),1,number_seeds);toannual(nominal_household_consumption_quarterly(5-q:end-mod(q,4),:)')'];
    model.nominal_household_consumption_growth=diff(log([repmat(data.nominal_household_consumption(data.years_num==year_num),1,number_seeds);toannual(nominal_household_consumption_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.nominal_household_consumption_growth=exp(model.nominal_household_consumption_growth)-1;
    model.nominal_household_consumption_growth=[repmat(data.nominal_household_consumption_growth(data.years_num==year_num),1,number_seeds);model.nominal_household_consumption_growth];
    
    model.nominal_household_consumption_growth_quarterly=diff(log([repmat(data.nominal_household_consumption_quarterly(data.quarters_num==quarter_num),1,number_seeds);nominal_household_consumption_quarterly]));
    % calculate discrete compounding rate
    model.nominal_household_consumption_growth_quarterly=exp(model.nominal_household_consumption_growth_quarterly)-1;
    model.nominal_household_consumption_growth_quarterly=[repmat(data.nominal_household_consumption_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.nominal_household_consumption_growth_quarterly];
    
    household_consumption_deflator_quarterly=nominal_household_consumption_quarterly./real_household_consumption_quarterly;
      
    model.household_consumption_deflator_quarterly=[repmat(data.household_consumption_deflator_quarterly(data.quarters_num==quarter_num),1,number_seeds);household_consumption_deflator_quarterly];
    model.household_consumption_deflator=[repmat(data.household_consumption_deflator(data.years_num==year_num),1,number_seeds);toannual_mean(household_consumption_deflator_quarterly(5-q:end-mod(q,4),:)')'];
    model.household_consumption_deflator_growth=diff(log([repmat(data.household_consumption_deflator(data.years_num==year_num),1,number_seeds);toannual_mean(household_consumption_deflator_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.household_consumption_deflator_growth=exp(model.household_consumption_deflator_growth)-1;
    model.household_consumption_deflator_growth=[repmat(data.household_consumption_deflator_growth(data.years_num==year_num),1,number_seeds);model.household_consumption_deflator_growth];
    
    model.household_consumption_deflator_growth_quarterly=diff(log([repmat(data.household_consumption_deflator_quarterly(data.quarters_num==quarter_num),1,number_seeds);household_consumption_deflator_quarterly]));
    % calculate discrete compounding rate
    model.household_consumption_deflator_growth_quarterly=exp(model.household_consumption_deflator_growth_quarterly)-1;
    model.household_consumption_deflator_growth_quarterly=[repmat(data.household_consumption_deflator_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.household_consumption_deflator_growth_quarterly];
    
    real_government_consumption_growth_quarterly=diff(log(real_government_consumption));
    real_government_consumption_growth_quarterly=exp(real_government_consumption_growth_quarterly)-1;
    real_government_consumption_quarterly=data.real_government_consumption_quarterly(data.quarters_num==quarter_num)*cumprod(1+real_government_consumption_growth_quarterly);
    
    model.real_government_consumption_quarterly=[repmat(data.real_government_consumption_quarterly(data.quarters_num==quarter_num),1,number_seeds);real_government_consumption_quarterly];
    model.real_government_consumption=[repmat(data.real_government_consumption(data.years_num==year_num),1,number_seeds);toannual(real_government_consumption_quarterly(5-q:end-mod(q,4),:)')'];
    model.real_government_consumption_growth=diff(log([repmat(data.real_government_consumption(data.years_num==year_num),1,number_seeds);toannual(real_government_consumption_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.real_government_consumption_growth=exp(model.real_government_consumption_growth)-1;
    model.real_government_consumption_growth=[repmat(data.real_government_consumption_growth(data.years_num==year_num),1,number_seeds);model.real_government_consumption_growth];
    
    model.real_government_consumption_growth_quarterly=diff(log([repmat(data.real_government_consumption_quarterly(data.quarters_num==quarter_num),1,number_seeds);real_government_consumption_quarterly]));
    % calculate discrete compounding rate
    model.real_government_consumption_growth_quarterly=exp(model.real_government_consumption_growth_quarterly)-1;
    model.real_government_consumption_growth_quarterly=[repmat(data.real_government_consumption_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.real_government_consumption_growth_quarterly];
    
    nominal_government_consumption_growth_quarterly=diff(log(nominal_government_consumption));
    nominal_government_consumption_growth_quarterly=exp(nominal_government_consumption_growth_quarterly)-1;
    nominal_government_consumption_quarterly=data.nominal_government_consumption_quarterly(data.quarters_num==quarter_num)*cumprod(1+nominal_government_consumption_growth_quarterly);
    
    model.nominal_government_consumption_quarterly=[repmat(data.nominal_government_consumption_quarterly(data.quarters_num==quarter_num),1,number_seeds);nominal_government_consumption_quarterly];
    model.nominal_government_consumption=[repmat(data.nominal_government_consumption(data.years_num==year_num),1,number_seeds);toannual(nominal_government_consumption_quarterly(5-q:end-mod(q,4),:)')'];
    model.nominal_government_consumption_growth=diff(log([repmat(data.nominal_government_consumption(data.years_num==year_num),1,number_seeds);toannual(nominal_government_consumption_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.nominal_government_consumption_growth=exp(model.nominal_government_consumption_growth)-1;
    model.nominal_government_consumption_growth=[repmat(data.nominal_government_consumption_growth(data.years_num==year_num),1,number_seeds);model.nominal_government_consumption_growth];
    
    model.nominal_government_consumption_growth_quarterly=diff(log([repmat(data.nominal_government_consumption_quarterly(data.quarters_num==quarter_num),1,number_seeds);nominal_government_consumption_quarterly]));
    % calculate discrete compounding rate
    model.nominal_government_consumption_growth_quarterly=exp(model.nominal_government_consumption_growth_quarterly)-1;
    model.nominal_government_consumption_growth_quarterly=[repmat(data.nominal_government_consumption_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.nominal_government_consumption_growth_quarterly];
    
    government_consumption_deflator_quarterly=nominal_government_consumption_quarterly./real_government_consumption_quarterly;
      
    model.government_consumption_deflator_quarterly=[repmat(data.government_consumption_deflator_quarterly(data.quarters_num==quarter_num),1,number_seeds);government_consumption_deflator_quarterly];
    model.government_consumption_deflator=[repmat(data.government_consumption_deflator(data.years_num==year_num),1,number_seeds);toannual_mean(government_consumption_deflator_quarterly(5-q:end-mod(q,4),:)')'];
    model.government_consumption_deflator_growth=diff(log([repmat(data.government_consumption_deflator(data.years_num==year_num),1,number_seeds);toannual_mean(government_consumption_deflator_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.government_consumption_deflator_growth=exp(model.government_consumption_deflator_growth)-1;
    model.government_consumption_deflator_growth=[repmat(data.government_consumption_deflator_growth(data.years_num==year_num),1,number_seeds);model.government_consumption_deflator_growth];
    
    model.government_consumption_deflator_growth_quarterly=diff(log([repmat(data.government_consumption_deflator_quarterly(data.quarters_num==quarter_num),1,number_seeds);government_consumption_deflator_quarterly]));
    % calculate discrete compounding rate
    model.government_consumption_deflator_growth_quarterly=exp(model.government_consumption_deflator_growth_quarterly)-1;
    model.government_consumption_deflator_growth_quarterly=[repmat(data.government_consumption_deflator_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.government_consumption_deflator_growth_quarterly];
    
    real_capitalformation_growth_quarterly=diff(log(real_capitalformation));
    real_capitalformation_growth_quarterly=exp(real_capitalformation_growth_quarterly)-1;
    real_capitalformation_quarterly=data.real_capitalformation_quarterly(data.quarters_num==quarter_num)*cumprod(1+real_capitalformation_growth_quarterly);
    
    model.real_capitalformation_quarterly=[repmat(data.real_capitalformation_quarterly(data.quarters_num==quarter_num),1,number_seeds);real_capitalformation_quarterly];
    model.real_capitalformation=[repmat(data.real_capitalformation(data.years_num==year_num),1,number_seeds);toannual(real_capitalformation_quarterly(5-q:end-mod(q,4),:)')'];
    model.real_capitalformation_growth=diff(log([repmat(data.real_capitalformation(data.years_num==year_num),1,number_seeds);toannual(real_capitalformation_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.real_capitalformation_growth=exp(model.real_capitalformation_growth)-1;
    model.real_capitalformation_growth=[repmat(data.real_capitalformation_growth(data.years_num==year_num),1,number_seeds);model.real_capitalformation_growth];
    
    model.real_capitalformation_growth_quarterly=diff(log([repmat(data.real_capitalformation_quarterly(data.quarters_num==quarter_num),1,number_seeds);real_capitalformation_quarterly]));
    % calculate discrete compounding rate
    model.real_capitalformation_growth_quarterly=exp(model.real_capitalformation_growth_quarterly)-1;
    model.real_capitalformation_growth_quarterly=[repmat(data.real_capitalformation_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.real_capitalformation_growth_quarterly];
    
    nominal_capitalformation_growth_quarterly=diff(log(nominal_capitalformation));
    nominal_capitalformation_growth_quarterly=exp(nominal_capitalformation_growth_quarterly)-1;
    nominal_capitalformation_quarterly=data.nominal_capitalformation_quarterly(data.quarters_num==quarter_num)*cumprod(1+nominal_capitalformation_growth_quarterly);
    
    model.nominal_capitalformation_quarterly=[repmat(data.nominal_capitalformation_quarterly(data.quarters_num==quarter_num),1,number_seeds);nominal_capitalformation_quarterly];
    model.nominal_capitalformation=[repmat(data.nominal_capitalformation(data.years_num==year_num),1,number_seeds);toannual(nominal_capitalformation_quarterly(5-q:end-mod(q,4),:)')'];
    model.nominal_capitalformation_growth=diff(log([repmat(data.nominal_capitalformation(data.years_num==year_num),1,number_seeds);toannual(nominal_capitalformation_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.nominal_capitalformation_growth=exp(model.nominal_capitalformation_growth)-1;
    model.nominal_capitalformation_growth=[repmat(data.nominal_capitalformation_growth(data.years_num==year_num),1,number_seeds);model.nominal_capitalformation_growth];
    
    model.nominal_capitalformation_growth_quarterly=diff(log([repmat(data.nominal_capitalformation_quarterly(data.quarters_num==quarter_num),1,number_seeds);nominal_capitalformation_quarterly]));
    % calculate discrete compounding rate
    model.nominal_capitalformation_growth_quarterly=exp(model.nominal_capitalformation_growth_quarterly)-1;
    model.nominal_capitalformation_growth_quarterly=[repmat(data.nominal_capitalformation_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.nominal_capitalformation_growth_quarterly];
    
    capitalformation_deflator_quarterly=nominal_capitalformation_quarterly./real_capitalformation_quarterly;
      
    model.capitalformation_deflator_quarterly=[repmat(data.capitalformation_deflator_quarterly(data.quarters_num==quarter_num),1,number_seeds);capitalformation_deflator_quarterly];
    model.capitalformation_deflator=[repmat(data.capitalformation_deflator(data.years_num==year_num),1,number_seeds);toannual_mean(capitalformation_deflator_quarterly(5-q:end-mod(q,4),:)')'];
    model.capitalformation_deflator_growth=diff(log([repmat(data.capitalformation_deflator(data.years_num==year_num),1,number_seeds);toannual_mean(capitalformation_deflator_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.capitalformation_deflator_growth=exp(model.capitalformation_deflator_growth)-1;
    model.capitalformation_deflator_growth=[repmat(data.capitalformation_deflator_growth(data.years_num==year_num),1,number_seeds);model.capitalformation_deflator_growth];
    
    model.capitalformation_deflator_growth_quarterly=diff(log([repmat(data.capitalformation_deflator_quarterly(data.quarters_num==quarter_num),1,number_seeds);capitalformation_deflator_quarterly]));
    % calculate discrete compounding rate
    model.capitalformation_deflator_growth_quarterly=exp(model.capitalformation_deflator_growth_quarterly)-1;
    model.capitalformation_deflator_growth_quarterly=[repmat(data.capitalformation_deflator_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.capitalformation_deflator_growth_quarterly];
    
    real_fixed_capitalformation_growth_quarterly=diff(log(real_fixed_capitalformation));
    real_fixed_capitalformation_growth_quarterly=exp(real_fixed_capitalformation_growth_quarterly)-1;
    real_fixed_capitalformation_quarterly=data.real_fixed_capitalformation_quarterly(data.quarters_num==quarter_num)*cumprod(1+real_fixed_capitalformation_growth_quarterly);
    
    model.real_fixed_capitalformation_quarterly=[repmat(data.real_fixed_capitalformation_quarterly(data.quarters_num==quarter_num),1,number_seeds);real_fixed_capitalformation_quarterly];
    model.real_fixed_capitalformation=[repmat(data.real_fixed_capitalformation(data.years_num==year_num),1,number_seeds);toannual(real_fixed_capitalformation_quarterly(5-q:end-mod(q,4),:)')'];
    model.real_fixed_capitalformation_growth=diff(log([repmat(data.real_fixed_capitalformation(data.years_num==year_num),1,number_seeds);toannual(real_fixed_capitalformation_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.real_fixed_capitalformation_growth=exp(model.real_fixed_capitalformation_growth)-1;
    model.real_fixed_capitalformation_growth=[repmat(data.real_fixed_capitalformation_growth(data.years_num==year_num),1,number_seeds);model.real_fixed_capitalformation_growth];
    
    model.real_fixed_capitalformation_growth_quarterly=diff(log([repmat(data.real_fixed_capitalformation_quarterly(data.quarters_num==quarter_num),1,number_seeds);real_fixed_capitalformation_quarterly]));
    % calculate discrete compounding rate
    model.real_fixed_capitalformation_growth_quarterly=exp(model.real_fixed_capitalformation_growth_quarterly)-1;
    model.real_fixed_capitalformation_growth_quarterly=[repmat(data.real_fixed_capitalformation_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.real_fixed_capitalformation_growth_quarterly];
    
    nominal_fixed_capitalformation_growth_quarterly=diff(log(nominal_fixed_capitalformation));
    nominal_fixed_capitalformation_growth_quarterly=exp(nominal_fixed_capitalformation_growth_quarterly)-1;
    nominal_fixed_capitalformation_quarterly=data.nominal_fixed_capitalformation_quarterly(data.quarters_num==quarter_num)*cumprod(1+nominal_fixed_capitalformation_growth_quarterly);
    
    model.nominal_fixed_capitalformation_quarterly=[repmat(data.nominal_fixed_capitalformation_quarterly(data.quarters_num==quarter_num),1,number_seeds);nominal_fixed_capitalformation_quarterly];
    model.nominal_fixed_capitalformation=[repmat(data.nominal_fixed_capitalformation(data.years_num==year_num),1,number_seeds);toannual(nominal_fixed_capitalformation_quarterly(5-q:end-mod(q,4),:)')'];
    model.nominal_fixed_capitalformation_growth=diff(log([repmat(data.nominal_fixed_capitalformation(data.years_num==year_num),1,number_seeds);toannual(nominal_fixed_capitalformation_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.nominal_fixed_capitalformation_growth=exp(model.nominal_fixed_capitalformation_growth)-1;
    model.nominal_fixed_capitalformation_growth=[repmat(data.nominal_fixed_capitalformation_growth(data.years_num==year_num),1,number_seeds);model.nominal_fixed_capitalformation_growth];
    
    model.nominal_fixed_capitalformation_growth_quarterly=diff(log([repmat(data.nominal_fixed_capitalformation_quarterly(data.quarters_num==quarter_num),1,number_seeds);nominal_fixed_capitalformation_quarterly]));
    % calculate discrete compounding rate
    model.nominal_fixed_capitalformation_growth_quarterly=exp(model.nominal_fixed_capitalformation_growth_quarterly)-1;
    model.nominal_fixed_capitalformation_growth_quarterly=[repmat(data.nominal_fixed_capitalformation_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.nominal_fixed_capitalformation_growth_quarterly];
    
    fixed_capitalformation_deflator_quarterly=nominal_fixed_capitalformation_quarterly./real_fixed_capitalformation_quarterly;
      
    model.fixed_capitalformation_deflator_quarterly=[repmat(data.fixed_capitalformation_deflator_quarterly(data.quarters_num==quarter_num),1,number_seeds);fixed_capitalformation_deflator_quarterly];
    model.fixed_capitalformation_deflator=[repmat(data.fixed_capitalformation_deflator(data.years_num==year_num),1,number_seeds);toannual_mean(fixed_capitalformation_deflator_quarterly(5-q:end-mod(q,4),:)')'];
    model.fixed_capitalformation_deflator_growth=diff(log([repmat(data.fixed_capitalformation_deflator(data.years_num==year_num),1,number_seeds);toannual_mean(fixed_capitalformation_deflator_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.fixed_capitalformation_deflator_growth=exp(model.fixed_capitalformation_deflator_growth)-1;
    model.fixed_capitalformation_deflator_growth=[repmat(data.fixed_capitalformation_deflator_growth(data.years_num==year_num),1,number_seeds);model.fixed_capitalformation_deflator_growth];
    
    model.fixed_capitalformation_deflator_growth_quarterly=diff(log([repmat(data.fixed_capitalformation_deflator_quarterly(data.quarters_num==quarter_num),1,number_seeds);fixed_capitalformation_deflator_quarterly]));
    % calculate discrete compounding rate
    model.fixed_capitalformation_deflator_growth_quarterly=exp(model.fixed_capitalformation_deflator_growth_quarterly)-1;
    model.fixed_capitalformation_deflator_growth_quarterly=[repmat(data.fixed_capitalformation_deflator_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.fixed_capitalformation_deflator_growth_quarterly];
    
    real_exports_growth_quarterly=diff(log(real_exports));
    real_exports_growth_quarterly=exp(real_exports_growth_quarterly)-1;
    real_exports_quarterly=data.real_exports_quarterly(data.quarters_num==quarter_num)*cumprod(1+real_exports_growth_quarterly);
    
    model.real_exports_quarterly=[repmat(data.real_exports_quarterly(data.quarters_num==quarter_num),1,number_seeds);real_exports_quarterly];
    model.real_exports=[repmat(data.real_exports(data.years_num==year_num),1,number_seeds);toannual(real_exports_quarterly(5-q:end-mod(q,4),:)')'];
    model.real_exports_growth=diff(log([repmat(data.real_exports(data.years_num==year_num),1,number_seeds);toannual(real_exports_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.real_exports_growth=exp(model.real_exports_growth)-1;
    model.real_exports_growth=[repmat(data.real_exports_growth(data.years_num==year_num),1,number_seeds);model.real_exports_growth];
    
    model.real_exports_growth_quarterly=diff(log([repmat(data.real_exports_quarterly(data.quarters_num==quarter_num),1,number_seeds);real_exports_quarterly]));
    % calculate discrete compounding rate
    model.real_exports_growth_quarterly=exp(model.real_exports_growth_quarterly)-1;
    model.real_exports_growth_quarterly=[repmat(data.real_exports_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.real_exports_growth_quarterly];
    
    nominal_exports_growth_quarterly=diff(log(nominal_exports));
    nominal_exports_growth_quarterly=exp(nominal_exports_growth_quarterly)-1;
    nominal_exports_quarterly=data.nominal_exports_quarterly(data.quarters_num==quarter_num)*cumprod(1+nominal_exports_growth_quarterly);
    
    model.nominal_exports_quarterly=[repmat(data.nominal_exports_quarterly(data.quarters_num==quarter_num),1,number_seeds);nominal_exports_quarterly];
    model.nominal_exports=[repmat(data.nominal_exports(data.years_num==year_num),1,number_seeds);toannual(nominal_exports_quarterly(5-q:end-mod(q,4),:)')'];
    model.nominal_exports_growth=diff(log([repmat(data.nominal_exports(data.years_num==year_num),1,number_seeds);toannual(nominal_exports_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.nominal_exports_growth=exp(model.nominal_exports_growth)-1;
    model.nominal_exports_growth=[repmat(data.nominal_exports_growth(data.years_num==year_num),1,number_seeds);model.nominal_exports_growth];
    
    model.nominal_exports_growth_quarterly=diff(log([repmat(data.nominal_exports_quarterly(data.quarters_num==quarter_num),1,number_seeds);nominal_exports_quarterly]));
    % calculate discrete compounding rate
    model.nominal_exports_growth_quarterly=exp(model.nominal_exports_growth_quarterly)-1;
    model.nominal_exports_growth_quarterly=[repmat(data.nominal_exports_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.nominal_exports_growth_quarterly];
    
    exports_deflator_quarterly=nominal_exports_quarterly./real_exports_quarterly;
      
    model.exports_deflator_quarterly=[repmat(data.exports_deflator_quarterly(data.quarters_num==quarter_num),1,number_seeds);exports_deflator_quarterly];
    model.exports_deflator=[repmat(data.exports_deflator(data.years_num==year_num),1,number_seeds);toannual_mean(exports_deflator_quarterly(5-q:end-mod(q,4),:)')'];
    model.exports_deflator_growth=diff(log([repmat(data.exports_deflator(data.years_num==year_num),1,number_seeds);toannual_mean(exports_deflator_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.exports_deflator_growth=exp(model.exports_deflator_growth)-1;
    model.exports_deflator_growth=[repmat(data.exports_deflator_growth(data.years_num==year_num),1,number_seeds);model.exports_deflator_growth];
    
    model.exports_deflator_growth_quarterly=diff(log([repmat(data.exports_deflator_quarterly(data.quarters_num==quarter_num),1,number_seeds);exports_deflator_quarterly]));
    % calculate discrete compounding rate
    model.exports_deflator_growth_quarterly=exp(model.exports_deflator_growth_quarterly)-1;
    model.exports_deflator_growth_quarterly=[repmat(data.exports_deflator_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.exports_deflator_growth_quarterly];
    
    real_imports_growth_quarterly=diff(log(real_imports));
    real_imports_growth_quarterly=exp(real_imports_growth_quarterly)-1;
    real_imports_quarterly=data.real_imports_quarterly(data.quarters_num==quarter_num)*cumprod(1+real_imports_growth_quarterly);
    
    model.real_imports_quarterly=[repmat(data.real_imports_quarterly(data.quarters_num==quarter_num),1,number_seeds);real_imports_quarterly];
    model.real_imports=[repmat(data.real_imports(data.years_num==year_num),1,number_seeds);toannual(real_imports_quarterly(5-q:end-mod(q,4),:)')'];
    model.real_imports_growth=diff(log([repmat(data.real_imports(data.years_num==year_num),1,number_seeds);toannual(real_imports_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.real_imports_growth=exp(model.real_imports_growth)-1;
    model.real_imports_growth=[repmat(data.real_imports_growth(data.years_num==year_num),1,number_seeds);model.real_imports_growth];
    
    model.real_imports_growth_quarterly=diff(log([repmat(data.real_imports_quarterly(data.quarters_num==quarter_num),1,number_seeds);real_imports_quarterly]));
    % calculate discrete compounding rate
    model.real_imports_growth_quarterly=exp(model.real_imports_growth_quarterly)-1;
    model.real_imports_growth_quarterly=[repmat(data.real_imports_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.real_imports_growth_quarterly];
    
    nominal_imports_growth_quarterly=diff(log(nominal_imports));
    nominal_imports_growth_quarterly=exp(nominal_imports_growth_quarterly)-1;
    nominal_imports_quarterly=data.nominal_imports_quarterly(data.quarters_num==quarter_num)*cumprod(1+nominal_imports_growth_quarterly);
    
    model.nominal_imports_quarterly=[repmat(data.nominal_imports_quarterly(data.quarters_num==quarter_num),1,number_seeds);nominal_imports_quarterly];
    model.nominal_imports=[repmat(data.nominal_imports(data.years_num==year_num),1,number_seeds);toannual(nominal_imports_quarterly(5-q:end-mod(q,4),:)')'];
    model.nominal_imports_growth=diff(log([repmat(data.nominal_imports(data.years_num==year_num),1,number_seeds);toannual(nominal_imports_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.nominal_imports_growth=exp(model.nominal_imports_growth)-1;
    model.nominal_imports_growth=[repmat(data.nominal_imports_growth(data.years_num==year_num),1,number_seeds);model.nominal_imports_growth];
    
    model.nominal_imports_growth_quarterly=diff(log([repmat(data.nominal_imports_quarterly(data.quarters_num==quarter_num),1,number_seeds);nominal_imports_quarterly]));
    % calculate discrete compounding rate
    model.nominal_imports_growth_quarterly=exp(model.nominal_imports_growth_quarterly)-1;
    model.nominal_imports_growth_quarterly=[repmat(data.nominal_imports_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.nominal_imports_growth_quarterly];
    
    imports_deflator_quarterly=nominal_imports_quarterly./real_imports_quarterly;
      
    model.imports_deflator_quarterly=[repmat(data.imports_deflator_quarterly(data.quarters_num==quarter_num),1,number_seeds);imports_deflator_quarterly];
    model.imports_deflator=[repmat(data.imports_deflator(data.years_num==year_num),1,number_seeds);toannual_mean(imports_deflator_quarterly(5-q:end-mod(q,4),:)')'];
    model.imports_deflator_growth=diff(log([repmat(data.imports_deflator(data.years_num==year_num),1,number_seeds);toannual_mean(imports_deflator_quarterly(5-q:end-mod(q,4),:)')']));
    % calculate discrete compounding rate
    model.imports_deflator_growth=exp(model.imports_deflator_growth)-1;
    model.imports_deflator_growth=[repmat(data.imports_deflator_growth(data.years_num==year_num),1,number_seeds);model.imports_deflator_growth];
    
    model.imports_deflator_growth_quarterly=diff(log([repmat(data.imports_deflator_quarterly(data.quarters_num==quarter_num),1,number_seeds);imports_deflator_quarterly]));
    % calculate discrete compounding rate
    model.imports_deflator_growth_quarterly=exp(model.imports_deflator_growth_quarterly)-1;
    model.imports_deflator_growth_quarterly=[repmat(data.imports_deflator_growth_quarterly(data.quarters_num==quarter_num),1,number_seeds);model.imports_deflator_growth_quarterly];
    
    operating_surplus_growth_quarterly=diff(log(operating_surplus));
    operating_surplus_growth_quarterly=exp(operating_surplus_growth_quarterly)-1;
    operating_surplus_quarterly=data.operating_surplus_quarterly(data.quarters_num==quarter_num)*cumprod(1+operating_surplus_growth_quarterly);
    
    model.operating_surplus_quarterly=[repmat(data.operating_surplus_quarterly(data.quarters_num==quarter_num),1,number_seeds);operating_surplus_quarterly];
    model.operating_surplus=[repmat(data.operating_surplus(data.years_num==year_num),1,number_seeds);toannual(operating_surplus_quarterly(5-q:end-mod(q,4),:)')'];
    
    compensation_employees_growth_quarterly=diff(log(compensation_employees));
    compensation_employees_growth_quarterly=exp(compensation_employees_growth_quarterly)-1;
    compensation_employees_quarterly=data.compensation_employees_quarterly(data.quarters_num==quarter_num)*cumprod(1+compensation_employees_growth_quarterly);
    
    model.compensation_employees_quarterly=[repmat(data.compensation_employees_quarterly(data.quarters_num==quarter_num),1,number_seeds);compensation_employees_quarterly];
    model.compensation_employees=[repmat(data.compensation_employees(data.years_num==year_num),1,number_seeds);toannual(compensation_employees_quarterly(5-q:end-mod(q,4),:)')'];
    
    wages_growth_quarterly=diff(log(wages));
    wages_growth_quarterly=exp(wages_growth_quarterly)-1;
    wages_quarterly=data.wages_quarterly(data.quarters_num==quarter_num)*cumprod(1+wages_growth_quarterly);
    
    model.wages_quarterly=[repmat(data.wages_quarterly(data.quarters_num==quarter_num),1,number_seeds);wages_quarterly];
    model.wages=[repmat(data.wages(data.years_num==year_num),1,number_seeds);toannual(wages_quarterly(5-q:end-mod(q,4),:)')'];
    
    model.quarters_num=datenum(datetime(year(datetime(datestr(quarter_num))),month(datetime(datestr(quarter_num)))+(0:3:3*horizon)+1,0))';
    model.years_num=datenum(year(datetime(datestr(quarter_num)))+1,1:12:horizon/4*12+floor(q/4),0)';
    
    model.nominal_nace10_gva_quarterly(:,:,1)=sum(nominal_sector_gva(:,:,1:3),3)*data.nominal_nace10_gva_quarterly(data.quarters_num==quarter_num,1)/sum(nominal_sector_gva(1,1,1:3));
    model.nominal_nace10_gva_quarterly(:,:,2)=sum(nominal_sector_gva(:,:,4:26),3)*data.nominal_nace10_gva_quarterly(data.quarters_num==quarter_num,2)/sum(nominal_sector_gva(1,1,4:26));
    model.nominal_nace10_gva_quarterly(:,:,3)=sum(nominal_sector_gva(:,:,27),3)*data.nominal_nace10_gva_quarterly(data.quarters_num==quarter_num,3)/sum(nominal_sector_gva(1,1,27));
    model.nominal_nace10_gva_quarterly(:,:,4)=sum(nominal_sector_gva(:,:,28:36),3)*data.nominal_nace10_gva_quarterly(data.quarters_num==quarter_num,4)/sum(nominal_sector_gva(1,1,28:36));
    model.nominal_nace10_gva_quarterly(:,:,5)=sum(nominal_sector_gva(:,:,37:40),3)*data.nominal_nace10_gva_quarterly(data.quarters_num==quarter_num,5)/sum(nominal_sector_gva(1,1,37:40));
    model.nominal_nace10_gva_quarterly(:,:,6)=sum(nominal_sector_gva(:,:,41:43),3)*data.nominal_nace10_gva_quarterly(data.quarters_num==quarter_num,6)/sum(nominal_sector_gva(1,1,41:43));
    model.nominal_nace10_gva_quarterly(:,:,7)=sum(nominal_sector_gva(:,:,44),3)*data.nominal_nace10_gva_quarterly(data.quarters_num==quarter_num,7)/sum(nominal_sector_gva(1,1,44));
    model.nominal_nace10_gva_quarterly(:,:,8)=sum(nominal_sector_gva(:,:,45:53),3)*data.nominal_nace10_gva_quarterly(data.quarters_num==quarter_num,8)/sum(nominal_sector_gva(1,1,45:53));
    model.nominal_nace10_gva_quarterly(:,:,9)=sum(nominal_sector_gva(:,:,54:57),3)*data.nominal_nace10_gva_quarterly(data.quarters_num==quarter_num,9)/sum(nominal_sector_gva(1,1,54:57));
    model.nominal_nace10_gva_quarterly(:,:,10)=sum(nominal_sector_gva(:,:,58:62),3)*data.nominal_nace10_gva_quarterly(data.quarters_num==quarter_num,10)/sum(nominal_sector_gva(1,1,58:62));
    
    nominal_nace10_gva_fints=toannual(fints(model.quarters_num(6-q:end-mod(q,4)),model.nominal_nace10_gva_quarterly(6-q:end-mod(q,4),:,:)),'CalcMethod','cumsum','BusDays',0);
    model.nominal_nace10_gva=[reshape(repmat(data.nominal_nace10_gva(data.years_num==model.years_num(1),:),[number_seeds,1]),1,number_seeds,10);reshape(fts2mat(nominal_nace10_gva_fints),length(model.years_num)-1,number_seeds,10)];
    
    model.real_nace10_gva_quarterly(:,:,1)=sum(real_sector_gva(:,:,1:3),3)*data.real_nace10_gva_quarterly(data.quarters_num==quarter_num,1)/sum(real_sector_gva(1,1,1:3));
    model.real_nace10_gva_quarterly(:,:,2)=sum(real_sector_gva(:,:,4:26),3)*data.real_nace10_gva_quarterly(data.quarters_num==quarter_num,2)/sum(real_sector_gva(1,1,4:26));
    model.real_nace10_gva_quarterly(:,:,3)=sum(real_sector_gva(:,:,27),3)*data.real_nace10_gva_quarterly(data.quarters_num==quarter_num,3)/sum(real_sector_gva(1,1,27));
    model.real_nace10_gva_quarterly(:,:,4)=sum(real_sector_gva(:,:,28:36),3)*data.real_nace10_gva_quarterly(data.quarters_num==quarter_num,4)/sum(real_sector_gva(1,1,28:36));
    model.real_nace10_gva_quarterly(:,:,5)=sum(real_sector_gva(:,:,37:40),3)*data.real_nace10_gva_quarterly(data.quarters_num==quarter_num,5)/sum(real_sector_gva(1,1,37:40));
    model.real_nace10_gva_quarterly(:,:,6)=sum(real_sector_gva(:,:,41:43),3)*data.real_nace10_gva_quarterly(data.quarters_num==quarter_num,6)/sum(real_sector_gva(1,1,41:43));
    model.real_nace10_gva_quarterly(:,:,7)=sum(real_sector_gva(:,:,44),3)*data.real_nace10_gva_quarterly(data.quarters_num==quarter_num,7)/sum(real_sector_gva(1,1,44));
    model.real_nace10_gva_quarterly(:,:,8)=sum(real_sector_gva(:,:,45:53),3)*data.real_nace10_gva_quarterly(data.quarters_num==quarter_num,8)/sum(real_sector_gva(1,1,45:53));
    model.real_nace10_gva_quarterly(:,:,9)=sum(real_sector_gva(:,:,54:57),3)*data.real_nace10_gva_quarterly(data.quarters_num==quarter_num,9)/sum(real_sector_gva(1,1,54:57));
    model.real_nace10_gva_quarterly(:,:,10)=sum(real_sector_gva(:,:,58:62),3)*data.real_nace10_gva_quarterly(data.quarters_num==quarter_num,10)/sum(real_sector_gva(1,1,58:62));
    
    real_nace10_gva_fints=toannual(fints(model.quarters_num(6-q:end-mod(q,4)),model.real_nace10_gva_quarterly(6-q:end-mod(q,4),:,:)),'CalcMethod','cumsum','BusDays',0);
    model.real_nace10_gva=[reshape(repmat(data.real_nace10_gva(data.years_num==model.years_num(1),:),[number_seeds,1]),1,number_seeds,10);reshape(fts2mat(real_nace10_gva_fints),length(model.years_num)-1,number_seeds,10)];
    
    model.euribor=(1+euribor).^4-1;
    
    save(['../data/abmx/',num2str(year(datetime(datestr(quarters_num(i))))),'Q',num2str(quarter(datetime(datestr(quarters_num(i))))),'.mat'],'model');    
end
