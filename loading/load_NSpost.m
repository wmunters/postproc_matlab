function [NSpost1] = load_NSpost(name)
NS = load(name);

NSpost1.t = NS(:,1);
NSpost1.Um = NS(:,2);
NSpost1.Vm = NS(:,3);
NSpost1.Wm = NS(:,4);
NSpost1.Etot = NS(:,5);
NSpost1.Eturb = NS(:,6);
NSpost1.frmx = NS(:,7);
NSpost1.frmy = NS(:,8);
NSpost1.frUm = NS(:,9);
if size(NSpost1,2)==10; 
    NSpost1.frsgs = NS(:,10);
end

