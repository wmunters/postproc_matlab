function BLtstatcc = load_BLtstatcc(filename)    
%Function to only import BLtstatcc, useful to compare different simulations
BL_tcc = load(filename);

fid = fopen(filename);
    A = fscanf(fid,'%s');
    [~,A] = strtok(A,'=');
    [n,A] = strtok(A(2:100),';');            %20 just arbitrary number 
    [~,A] = strtok(A,'=');
    [tau_av,~] = strtok(A(2:20),'%');
    n = str2num(n);
    tau_av = str2num(tau_av);
    
    BLtstatcc.n  = n;                    %number of samples
    BLtstatcc.tau_av = tau_av;           %average shear stress
    BLtstatcc.z  = BL_tcc(:,1);
    BLtstatcc.Um = BL_tcc(:,2);
    BLtstatcc.Vm = BL_tcc(:,3);
    BLtstatcc.Wm = BL_tcc(:,4);
    BLtstatcc.uu = BL_tcc(:,5);
    BLtstatcc.vv = BL_tcc(:,6);
    BLtstatcc.ww = BL_tcc(:,7);
    BLtstatcc.uv = BL_tcc(:,8);
    BLtstatcc.uw = BL_tcc(:,9);
    BLtstatcc.vw = BL_tcc(:,10);
end