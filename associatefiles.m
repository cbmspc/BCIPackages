commandwindow; 
cwd=pwd; 
cd([matlabroot '\toolbox\matlab\winfun\private']); 
fileassoc('add',{'.m','.mat','.fig','.p','.mdl',['.' mexext]}); 
cd(cwd);
