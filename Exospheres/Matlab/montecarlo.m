% Event-driver for ballistically hopping molecules
% each particle has:
%    position p_r (longitude,latitude),
%    status (on surface or in flight),
%    time (until it arrives on surface or until it will leave the surface)


function [p_r,p_s,p_t,p_n,ccc] = montecarlo(NP,p_r,p_s,p_t,p_n,Tsurf,dtsec,ccc,Qdissoc)
  % particle migration within a time interval of duration dtsec
  % event-driven 
  
  for i=1:NP
    % sprintf('montecarlo: working on particle %d',i)

    while p_t(i) <= dtsec  % do this for a duration of dtsec
      if p_s(i)<0,
	break  % exit while loop and move on to next particle
      end
      
      % sprintf('next event is %d %f %f',i,p_t(i),p_s(i))
      switch p_s(i)
        case num2cell(-10:-1)
          break
        case{0}  % leaving surface
          [ii,jj] = inbox(p_r(i,:));
          [p_r(i,:), p_s(i), p_t(i)] = hop1(p_r(i,:),p_t(i),Tsurf(ii,jj),Qdissoc);
          if p_s(i)==-1, ccc(1)=ccc(1)+1; end
          if p_s(i)==-2, ccc(2)=ccc(2)+1; end
          p_n(i) = p_n(i)+1;  % count hops
        case{1}  % landing
          [ii,jj] = inbox(p_r(i,:));
          if incoldtrap(p_r(i,:)),
            if p_r(i,2)>0., p_s(i)=-3; end
            if p_r(i,2)<0., p_s(i)=-4; end
            if p_s(i)==-3, ccc(3)=ccc(3)+1; end
            if p_s(i)==-4, ccc(4)=ccc(4)+1; end
            p_t(i) = residence_time(90.);  % very long
	  else
            residencetime = residence_time(Tsurf(ii,jj));
            p_t(i) = p_t(i) + residencetime;
            p_s(i) = 0;
	  end
	  %disp(sprintf("M %d %f %f %d %g %d\n",i,p_r(i,:),p_s(i),p_t(i),p_n(i)))
      end
    end   % end of time/event loop

  end  % end loop over particles

  % subtract dtsec from all times
  k = find( p_s>=0 );
  p_t(k) = p_t(k) - dtsec;
  
end 
