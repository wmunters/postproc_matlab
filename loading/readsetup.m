function setup = readsetup

k = 1;
fid = fopen('NS.setup');
while ~feof(fid)
  results(k,:) = textscan(fid, '%s %s ', 1,'headerlines',1);
  fgets(fid);
  k = k+1;
end
% Has to be adapted in case of high Re wall model 
setup.casename = results{1,1}{1};
setup.casename = setup.casename(1:length(setup.casename)-1);
setup.Nx = str2double(results{2,1}{1});
setup.Ny = str2double(results{2,2}{1});
setup.Lx = str2double(strtok(results{3,1}{1},'d')); % CAN ONLY BE USED IF EXPONENT AFTER d IS ZERO
setup.Ly = str2double(strtok(results{3,2}{1},'d'));
setup.Re = str2double(strtok(results{4,1}{1},'d'));
setup.LRe = str2double(strtok(results{5,1}{1},'d'));
setup.URe = str2double(strtok(results{6,1}{1},'d'));
setup.tstart = results{7,1}{1};
setup.tstop = results{8,1}{1};
setup.tpost1 = results{9,1}{1};
setup.tpost2 = results{10,1}{1};
setup.tpost3 = results{11,1}{1};
setup.meshname = results{12,1}{1};
setup.boundaryname = results{13,1}{1};
%   setup.boundaryname = setup.boundaryname(1:length(setup.boundaryname)-1); % Weet niet waarom dit juist moet. Vreemd
setup.IBtype = results{14,1}{1};
%setup.initialcondname = results(14,1){1};
setup.bctype = str2double(results{16,1}{1});
setup.model = str2double(results{17,1}{1});
setup.smagcoeff = results{18,1}{1};
setup.dealias = results{19,1}{1};
setup.CFLconv = results{20,1}{1};
setup.CFLvisc = results{21,1}{1};
setup.stat = results{22,1}{1};

end