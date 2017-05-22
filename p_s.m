function [z] = p_s(dlist,time,start)
% Get prompts per second/100
time = ceil(time*0.1);
p=0;
T=start;
tmplen=0;
z = zeros(time,1);
for i=1:length(dlist)
  if (dlist(i)<(2^31))&&(dlist(i)>=(2^30))
    p=p+1;
  elseif (dlist(i)>=(2^31))&&(dlist(i)<2684354560)
    if mod(T,10)==0
      tmptime = T*0.1;
      if tmptime>0
        z(tmptime) = p-tmplen;
        tmplen=p;
      end
    end
    T=T+1;
  end
end

end

