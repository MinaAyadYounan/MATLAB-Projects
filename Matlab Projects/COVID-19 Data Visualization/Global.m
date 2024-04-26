classdef Global 
properties
countries    % there is class called countries for every country and Global saved in this variable 
countriesName % name of all countries used to add countries as items in list box
              % startDate and EndDate for making the date dynamic   
startDate      
endDate
end
properties(Access = private)
comulativeInjured % used this variable to sum all comulative injured in each country add pass it to class country with country name = 'Global'; 
comulativeDead    % sum all comultive deaths
     % note : all calculations on comulativeInjured and comulativeDead in States class 
end
methods
    function obj = Global()
     load('iKxt-dBOQfisbfnQTkH48A_37ac1662db254c17841c5c0cc420eef1_covid_data.mat','covid_data');
     allcountries = covid_data(2:end,1);                        % get all countries    
     states = covid_data(2:end,2);                              % save all stastes in cloumn vector called states
     numericalArr = covid_data(2:end,3:end);                    % get all the numerical data for calcultions (comulative deaths and cases)
     unique_country = unique(covid_data(2:end,1));              % create variable have all country without any repetition  
     obj.startDate = covid_data(1,3);
     obj.endDate = covid_data(1,end);
     obj.countriesName =[ 'Global';unique_country];
     obj.countries{1} = 'Global';
     numericalArr= cell2mat(numericalArr); % covert cellArray to Array
     injured = numericalArr(:,1:2:end);    % injured cases have odd indices 
     dead = numericalArr(:,2:2:end);
     obj.comulativeInjured = zeros(1,length(covid_data(1,3:end))); % preallocation
     obj.comulativeDead =  zeros(1,length(covid_data(1,3:end)));
      % for loop for insert countries in country class
       for i = 1:length(unique_country)   
           index = count(allcountries,unique_country{i}); % to get the index of row of all countries have the same name to get the hole data in this country
           statesOfThiscountry = states(logical(index));  % stored states of this country
           injuredInCountry = injured(logical(index),:);
           injuredInDead = dead(logical(index),:);
             obj.comulativeInjured =obj.comulativeInjured +injuredInCountry(1,:) ; % sum all injured cases in each country 
             obj.comulativeDead =obj.comulativeDead + injuredInDead(1,:);
           
       % passing name of country , all states ,comulatives injured & dead in the conutry and it'as states 
            obj.countries{i+1} = Country(unique_country{i},statesOfThiscountry,injuredInCountry,injuredInDead);
       end
       obj.countries{1} = Country('Global',{'All'},obj.comulativeInjured,obj.comulativeDead);
    end
end
end