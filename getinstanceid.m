function instanceid = getinstanceid ()
% Instance ID = Computer ID + Process ID
cid = getcomputerid();
pid = getpid();
instanceid = [dec2base(cid, 36, 10) dec2base(pid, 36, 7)];

