function t = unixtime ()
utc_time = java.lang.System.currentTimeMillis;
t = utc_time / 1000;
