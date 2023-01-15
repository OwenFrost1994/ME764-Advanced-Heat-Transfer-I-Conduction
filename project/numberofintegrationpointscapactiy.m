%====================== No. integration points  for mass matrices =====================
%

function n = numberofintegrationpointscapactiy(nDof,nNoEl)
  
   if (nDof == 1)
       n = nNoEl;
   else if (nDof == 2)
           if (nNoEl == 3)
               n = 4;
           end
           if (nNoEl == 6)
               n = 7;
           end
           if (nNoEl == 4)
               n = 4;
           end
           if (nNoEl == 8)
               n = 9;
           end
       else if (nDof == 3)
               if (nNoEl == 4)
                   n = 4 ;
               end
               if (nNoEl == 10)
                   n = 5;
               end
               if (nNoEl == 8)
                   n = 27;
               end
               if (nNoEl == 20)
                   n = 27;
               end
           end
       end
   end
end