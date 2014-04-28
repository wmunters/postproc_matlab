function [BLst] = load_stat(name,Nx,Ny)

% First load in data
	BLd = load(name);
	BLd = reshape(BLd',[size(BLd,1)*size(BLd,2) 1]);
	
% Start dividing and reshaping, put into structure
    Nz = 80;
	N = Nx*Ny*Nz;
	
	BLst.u = reshape(BLd(0*N+1:N),[Nx Ny Nz]);
	BLst.v = reshape(BLd(1*N+1:2*N),[Nx Ny Nz]);
	BLst.w = reshape(BLd(2*N+1:3*N),[Nx Ny Nz]);
	BLst.uu = reshape(BLd(3*N+1:4*N),[Nx Ny Nz]);
	BLst.vv = reshape(BLd(4*N+1:5*N),[Nx Ny Nz]);
	BLst.ww = reshape(BLd(5*N+1:6*N),[Nx Ny Nz]);
	BLst.uv = reshape(BLd(6*N+1:7*N),[Nx Ny Nz]);
	BLst.uw = reshape(BLd(7*N+1:8*N),[Nx Ny Nz]);
	BLst.vw = reshape(BLd(8*N+1:9*N),[Nx Ny Nz]);
	BLst.p = reshape(BLd(9*N+1:10*N),[Nx Ny Nz]);
	BLst.w_st = reshape(BLd(10*N+1:11*N),[Nx Ny Nz]);


	
	

end
	