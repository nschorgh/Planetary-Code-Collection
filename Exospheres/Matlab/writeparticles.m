function writeparticles(unit,NP,p_r,p_s,p_t,p_n)
  for i=1:NP
    %[ii,jj] = inbox(p_r(i,:));
    fprintf(unit,"%d %f %f %d %g %d\n",i,p_r(i,:),p_s(i),p_t(i),p_n(i));
  end
end 
