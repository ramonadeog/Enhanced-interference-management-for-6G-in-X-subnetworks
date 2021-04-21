function [keeperX] = randmin(Npoints,mindist)
x = rand(1, 100000);
% Initialize first point.
keeperX = x(1);
% Try dropping down more points.
counter = 2;
k = 2;
while counter <= Npoints
  % Get a trial point.
  thisX = x(k);
  % See how far is is away from existing keeper points.
  distances = sqrt((thisX-keeperX).^2);
  minDistance = min(distances);
  if minDistance >= mindist
    keeperX(counter) = thisX;
    counter = counter + 1;
  end
  k = k+1;
end

end