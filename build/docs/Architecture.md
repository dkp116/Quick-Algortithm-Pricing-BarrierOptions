How the code will be layed


IDynamics  - picks the dyanmics that will be used, this will take in risk free rate and volatility 
eg BSM, Merton Jump Diffussion

Stock  -  Construced using the dynamics and a start price

IPricing - prices the derivative depending on the option and stock that is inputted and iterations is number of times we want montecarlo to run

IPricing will have the different pricing methods such as taylorapprox, closed form etc 

For certain combinations it will not work like taylor approx and bsm for that we will return null.


How we will implement these things eg Taylorapprox(stock (bsm) , barrier ) this should be invalid but 
Taylorapprox(stock (mjd) , barrier) should have an implementation, what check can we give we can use a factory method, when creating this centralises the process

So implmentation of IPricing will be Taylorapprox(stock (mjd) , barrier)  . price   -- add in the pricing logic stock.simulate this would use either bb or 1 increment 
then we can add in certain logical/ densities. i do not think i would use the 1 time increment for mjd i see no benefit to it. What if i want to use it however

For this purpose i think just use brownianbridge function and simulate function anad make the bb function null, gets the job done.

