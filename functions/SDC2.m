function [out1,out2]=SDC2(inp1, inp2, iconv, LatOrig, LonOrig)
% [out1,out2]=SDC2(inp1, inp2, iconv)
% Conversion between coordinate systems based on defined origin in (lat,lon) and rotation angle
% using short distance conversion.
% if iconv = 1: Cartesian to geographic:
%    inp1 = x, inp2 = y
%    out1 = lat, out2 = lon
% if iconv = -1: Geographic to cartesian:
%    inp1 = lat, inp2 = lon
%    out1 = x, out2 = y
%
global rearth ellip rlatc rad olat olon aa bb bc sint cost rotat 
% global LatOrig LonOrig rota

ifil=0; % set > 0 for verbose

% LatOrig = 50.2116;
% LonOrig = -126.5971;
% rota=-35;
% LatOrig = 54.47918;
% LonOrig = -84.91262;
rota=0;

orlat=LatOrig; % 50.2116;
orlon=LonOrig; %-126.5971;
% rota=-35;




setorg(orlat, orlon, rota, ifil);


if (iconv==1) % cartesian to geographic
    x = inp1; y = inp2;
    [xlat,xlon]=redist(x,y);
    out1=xlat; out2=xlon;
elseif (iconv==-1) % geographic to cartesian
    xlat = inp1; xlon = inp2;
    [x,y]=dist(xlat,xlon);
    out1=x; out2=y;
else
    error('Wrong value for iconv')
end

end

function [ xkm, ykm] = dist(xlat, xlon)
% Convert latitude and longitude to kilometers relative
% to center of coordinates by short distance conversion.
global rlatc rad olat olon aa bb sint cost rotat

% Set up short distance conversion by subr. SETORG
q = 60 .* xlat - olat;
yp = q + olat;
lat1 = atan(rlatc .* tan(rad .* yp ./ 60.0));
lat2 = atan(rlatc .* tan(rad .* olat ./ 60.0));
lat3 = (lat2 + lat1) ./ 2;
xx = 60 .* xlon - olon;  %  - wegen LON E
q = q .* aa;
xx = xx .* bb .* cos(lat3);
if rotat ~= 0
    % rotate coordinate system anticlockwise
    yp = cost .* q + sint .* xx;
    xx = cost .* xx - sint .* q;
    q = yp;
end

xkm = xx;
ykm = q;

end

function [ xlat, xlon] = redist(xkm, ykm)
% Convert from local Cartesian coordinates to latitude & longitude
global rlatc rad olat olon aa bb sint cost

xx = xkm;
yy = ykm;

% Rotate coordinates anticlockwise back
y = yy .* cost - xx .* sint;
x = yy .* sint + xx .* cost;
if (abs(aa) < 0.0000001)
    disp([' subr. redist: aa=',num2str(aa),', bb=',num2str(bb),', cos(lat1)=',num2str(clat1),'. Division by zero, run stops here'])
    error('redist>>> division by zero!!')
end
q = y ./ aa;
lat = (q+olat) ./ 60;
xlat = q+olat - 60 .* lat;
yp = 60 .* lat+xlat;
lat1 = atan(rlatc .* tan(yp .* rad ./ 60.0));
lat2 = atan(rlatc .* tan(olat .* rad ./ 60.));
lat3 = (lat1+lat2) ./ 2;
clat1 = cos(lat3);
bcl = bb .* clat1;
if (abs(bcl) < 0.000001)
    disp([' subr. redist: aa=',num2str(aa),', bb=',num2str(bb),', cos(lat1)=',num2str(clat1),'. Division by zero, run stops here'])
    error('redist>>> division by zero!!')
end
p = x ./ (bb .* clat1);
lon = (p+olon) ./ 60;
xlon = p+olon - 60 .* lon;
xlat = lat+xlat ./ 60;
xlon = lon+xlon ./ 60;

end


function setorg(orlat, orlon, rota, ifil)

global rearth ellip rlatc rad olat olon aa bb bc sint cost rotat

% O(r)LAT & O(r)LON : origin of cartesian coordinate system north west
%  if orlat and orlon are both set to zero , the Swiss Cartesian
%  coordinate system will be used (this system cannot be rotated).
%  For all other cases, orlat and orlon denote the origin of
%  the SDC.

rad = 0.017453292;

if (orlat==0.0 && orlon==0.0)
    olat = 46.95240;  % BERN (Switzerland) North
    olon = -7.439583 ; % BERN West
    rotat = 0.0;
else
    olat = orlat;
    olon = orlon;
    rotat = rota*rad;
end

olat = olat*60; % minutes N
olon = olon*60; % minutes W

%  NEW ELLIPSOID FOR WHOLE EARTH:   WGS72 == WORLD GEODETIC SYSTEM 1972

%  ALSO SET RLATC ACCORDING TO ORIGIN

rearth = 6378.135;
ellip = 298.26 ;        % (flattening)


%  CALCULATE RLATC:  CONVERSION FROM GEOGRAPHICAL LAT TO GEOCENTRICAL LAT

phi = olat*rad/60.0   ;           %  phi = geogr. lat
beta = phi-sin(phi*2.)/ellip;
rlatc = tan(beta)/tan(phi);

%  WRITE ELLIPSOIDE CONSTANTS

if(ifil>0)
    disp(['SHORT DISTANCE CONVERSION on ELLIPSOIDE of',...
        ' WORLD GEODETIC SYSTEM 1972 (WGS72)',...
        '==========================================',...
        '===================================',...
        '( Radius at equator (REARTH)= ',num2str(rearth),'km',...
        '(   1. / (ellipticity)      = ', num2str( ellip),...
        'Origin of cartesian coordinates [degrees]:'])
    if (orlat==0.0 && orlon==0.0)
        disp([' SWISS COORDINATE SYSTEM (we have to be special)',...
            ' (Origin = city of BERNE, Switzerland),  no rotation of grid, pos. y-axis toward N, pos. x-axis toward E'])
        
    else
        disp(['Origin: ',num2str(olat/60),' N, ',num2str(olon/60),' W '])
        disp([' without rotation of grid, pos. x-axis toward WEST, pos. y-axis toward NORTH'])
        disp(['Rotation of y-axis from North anticlock-wise with pos. angle given in degrees'])
        if rota > 0
            disp([' Rotation of grid anticlockwise by', num2str(rota),' degrees'])
        else
            disp([' Rotation of grid clockwise by',num2str(-1.*rota),' degrees'])
        end
    end
end

%   calculate aa &  bb
%   length of one minute of lat and lon in km

lat1 = atan(rlatc*tan(olat*rad/60.0)) ;      % geoc. lat for OLAT
lat2 = atan(rlatc*tan((olat+1.)*rad/60.0)) ; % geoc. lat for (OLAT+1min)
dela = lat2 - lat1;
r = rearth*(1.0 - (sin(lat1)^2)/ellip);      % kugelradius fuer lat=OLAT
aa = r*dela;   % aa = 1 min geogr. lat
delb = acos(sin(lat1)^2 + cos(rad/60.)*cos(lat1)^2);
bc = r*delb ;    % bc = 1 min geogr. lon
bb = r*delb/cos(lat1);

if(ifil>0)
    disp(['( Radius of sphere at OLAT = ',num2str(r),' km)'])
    disp(['Conversion of GEOGRAPHICAL LATITUDE to GEOCENTRICAL LATITUDE:'])
    disp(['RLATC = tan(GEOCENTR.LAT) / TAN(GEOGRAPH.LAT)','  RLATC = ',num2str(rlatc)])
    disp('Short distance conversions:')
    disp(['one min lat ~ ', num2str(aa),' km '])
    disp(['one min lon ~ ', num2str(bb),' km '])
    
end

%convert coordinates with rotation cosines (stored in Common)

sint = sin(rotat);
cost = cos(rotat);


end



