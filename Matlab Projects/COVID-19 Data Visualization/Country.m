classdef Country
 properties
 name    % store name of the country
 statesOfCountry % store states in class called States
 statesName % store all states of the country
 end
 methods
     function node = Country(countryName,states,injured,dead)
      node.name = countryName;
      len = size(states,1);
      states{1} = 'All'; % allocate the first element all and injured and dead belongs to the country 
            for i = 1:len
                node.statesName{i} = states{i};
       % passing name of state , injured, dead
            node.statesOfCountry{i} = States(states{i},injured(i,:),dead(i,:));
            end     
     end
 end
end