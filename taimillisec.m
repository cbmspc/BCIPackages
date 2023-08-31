function t = taimillisec ()
t = mod((unixtime+378691200)*1000,60000);
%  mod(unixtime + 378691200,60)*1000;