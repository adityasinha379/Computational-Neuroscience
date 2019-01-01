function d=VPSDM(tli,tlj,q)
nspi=length(tli);
nspj=length(tlj);

if q==0
   d=abs(nspi-nspj);
   return
elseif q==Inf
   d=nspi+nspj;
   return
end

scr=zeros(nspi+1,nspj+1);
scr(:,1)=(0:nspi)';
scr(1,:)=(0:nspj);
if(nspi && nspj)
   for i=2:nspi+1
      for j=2:nspj+1
         scr(i,j)=min([scr(i-1,j)+1 scr(i,j-1)+1 scr(i-1,j-1)+q*abs(tli(i-1)-tlj(j-1))]);
      end
   end
end
d=scr(nspi+1,nspj+1);
end