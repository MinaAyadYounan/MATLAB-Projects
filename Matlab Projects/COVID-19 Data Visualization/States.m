classdef States 
    properties
    st_name % states name 
    comulativeInjured 
    dailyInjured 
    comulativeDead
    dailyDead
    end
methods
    function this = States(name,inj,dead)
        this.st_name =name;
        this.comulativeInjured= inj(1,:);
        this.comulativeDead = dead(1,:);
        %  --------- calcultions ------------
        this.dailyInjured = [this.comulativeInjured(1),diff(this.comulativeInjured)]; % diff is built in function  subtract(vector(i+1) - vector(i))
        this.dailyInjured(this.dailyInjured < 0) = 0;  % assign all negative number  to Zero
        this.dailyDead = [this.comulativeDead(1),diff(this.comulativeDead)];
        this.dailyDead(this.dailyDead < 0) = 0; 
    end

   end
end