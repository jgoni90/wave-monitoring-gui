function [cp,cpp,freq,theta] = getParamValuesFromCP(data)

%Polar coordinates
[theta,rho] = cart2pol(data.arrowdata.p0(1),data.arrowdata.p0(2));

%Point correction
if rho > pi, rho = pi; end
if theta < 0
    theta = data.arrowdata.prev_p0p(1);
    rho = data.arrowdata.prev_p0p(2);
elseif theta > data.arrowdata.thetalim1 
    theta = data.arrowdata.thetalim1;
    rho = data.arrowdata.prev_p0p(2);
elseif theta < data.arrowdata.thetalim0
    theta = data.arrowdata.thetalim0;
    rho = data.arrowdata.prev_p0p(2);
end

%Save point correction
cpp = [theta,rho];
[cp(1),cp(2)] = pol2cart(theta,rho);

%Frequency map
freq = rho*data.arrowdata.a0 + data.arrowdata.b0;