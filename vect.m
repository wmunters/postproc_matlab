function [v] = vect(x)

   v= reshape(x',[size(x,1)*size(x,2) 1]);
   