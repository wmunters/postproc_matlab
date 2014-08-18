function field = load_field(filename,Nx,Ny,Nz,Nl)

   field = load(filename);
   field = reshape(field',[numel(field) 1]);
   if(nargin>4)
        field = reshape(field,[Nx Ny Nz Nl]);
   else
        field = reshape(field,[Nx Ny Nz]);
   end
   
