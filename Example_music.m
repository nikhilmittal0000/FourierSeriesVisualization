clc; clear; close all;

dt = 0.0005;
n = 5;
data_points = [-10.9575569358178,-50.9161490683232;
-10.2329192546588,-40.7712215320912;
-10.2329192546588,-29.5393374741202;
-10.2329192546588,-18.6697722567282;
-10.2329192546588,-3.81469979296017;
-10.5952380952378,10.6780538302279;
-9.87060041407885,24.4461697722565;
-9.87060041407885,40.3881987577636;
-10.2329192546588,54.880952380952;
-10.2329192546588,67.5621118012419;
-10.2329192546588,82.41718426501;
-9.50828157349885,96.185300207039;
-10.2329192546588,107.05486542443;
-10.2329192546588,122.272256728778;
-3.34886128364384,122.272256728778;
-1.53726708074585,107.77950310559;
1.72360248447215,97.6345755693578;
6.07142857142814,88.5766045548651;
12.9554865424432,79.5186335403723;
26.3612836438922,64.6635610766042;
38.3178053830232,48.3592132505173;
44.4772256728782,38.2142857142855;
47.3757763975152,29.1563146997926;
48.8250517598342,20.4606625258798;
50.2743271221532,13.5766045548658;
50.2743271221532,6.69254658385084;
48.4627329192542,0.533126293995835;
46.2888198757762,-8.88716356107618;
43.0279503105591,-19.0320910973082;
39.7670807453411,-28.4523809523812;
35.0569358178052,-38.5973084886132;
31.4337474120081,-44.7567287784682;
32.8830227743272,-37.1480331262942;
38.3178053830232,-25.1915113871632;
41.9409937888202,-12.5103519668732;
44.1149068322982,0.170807453415833;
44.1149068322982,11.7650103519668;
40.8540372670812,24.4461697722565;
37.5931677018632,33.1418219461694;
31.0714285714281,42.1997929606622;
20.9265010351962,50.5331262939955;
12.9554865424432,54.5186335403723;
7.52070393374716,57.7795031055897;
2.08592132505215,60.3157349896477;
-2.26190476190484,61.4026915113868;
-1.89958592132484,57.0548654244303;
-1.89958592132484,49.8084886128361;
-1.53726708074585,29.5186335403723;
-1.89958592132484,9.95341614906881;
-2.98654244306384,-21.5683229813662;
-2.62422360248485,-62.5103519668732;
-5.88509316770184,-74.8291925465842;
-10.9575569358178,-83.5248447204972;
-16.7546583850928,-89.6842650103522;
-24.7256728778464,-94.7567287784682;
-31.9720496894406,-98.7422360248442;
-42.1169772256725,-100.916149068323;
-51.8995859213247,-100.553830227743;
-59.8706004140783,-98.0175983436852;
-66.3923395445131,-92.5828157349892;
-69.6532091097305,-83.8871635610762;
-68.928571428571,-74.8291925465842;
-64.5807453416145,-67.2204968944102;
-60.232919254658,-61.7857142857142;
-51.5372670807449,-55.6262939958592;
-46.1024844720493,-52.3654244306422;
-34.8706004140783,-47.6552795031052;
-26.8995859213246,-47.6552795031052;
-17.8416149068328,-48.0175983436852;
-10.9575569358178,-50.9161490683232];

%Calculate all the fourier series coefficients
C0 = zeros(1,2);
C_positive = zeros(n,2);
C_negative = zeros(n,2);
[C0(1,1), C0(1,2)] = calc_Cn(data_points,0);
for i=1:n
    [C_positive(i,1), C_positive(i,2)] = calc_Cn(data_points, i);
    [C_negative(i,1), C_negative(i,2)] = calc_Cn(data_points, -i);
end

% Calculate function Values at Every time step
noOfPoints = 1/dt;
f = zeros(noOfPoints,2);
T = linspace(0,1,noOfPoints);
for i = 1:noOfPoints
    t = T(i);
    f_real = C0(1);
    f_imaginary = C0(2);
    for j = 1:n
        f_real = f_real + getFourierRealTerm(C_positive(j,:),j,t);
        f_imaginary = f_imaginary + getFourierImagTerm(C_positive(j,:),j,t);
    end
    for j = 1:n
        f_real = f_real + getFourierRealTerm(C_negative(j,:),-j,t);
        f_imaginary = f_imaginary + getFourierImagTerm(C_negative(j,:),-j,t);
    end
    f(i,:) = [f_real, f_imaginary];
end

% Calculate components for every vector/ arrow
Components0 = zeros(1,2,noOfPoints);
Components_positive = zeros(n,2,noOfPoints);
Components_negative = zeros(n,2,noOfPoints);
Vector_position0 = zeros(1,4,noOfPoints);
Vector_position_positive = zeros(n,4,noOfPoints);
Vector_position_negative = zeros(n,4,noOfPoints);

for k=1:noOfPoints
    t = T(k);
    Components0(1,:,k) = [C0(1) C0(2)];
    Vector_position0(1,:,k) = [0 0 C0(1) C0(2)];
    lastVectorPosition = [C0(1) C0(2)];
    for j = 1:n
        Components_positive(j,:,k) = [getFourierRealTerm(C_positive(j,:),j,t), getFourierImagTerm(C_positive(j,:),j,t)];
        Vector_position_positive(j,:,k) = [lastVectorPosition(1), lastVectorPosition(2), Components_positive(j,1,k)+lastVectorPosition(1), Components_positive(j,2,k)+lastVectorPosition(2)];
        lastVectorPosition = [Components_positive(j,1,k)+lastVectorPosition(1), Components_positive(j,2,k)+lastVectorPosition(2)];
        
        Components_negative(j,:,k) = [getFourierRealTerm(C_negative(j,:),-j,t), getFourierImagTerm(C_negative(j,:),-j,t)];
        Vector_position_negative(j,:,k) = [lastVectorPosition(1), lastVectorPosition(2), Components_negative(j,1,k)+lastVectorPosition(1), Components_negative(j,2,k)+lastVectorPosition(2)];
        lastVectorPosition = [Components_negative(j,1,k)+lastVectorPosition(1), Components_negative(j,2,k)+lastVectorPosition(2)];
    end
end

%Plot function values
functionH = plot(f(1,1),f(1,2),'Color',[1,1,0]);
xlim([-150 150]);ylim([-150 150]);
set(gca,'Color','k')
axis equal;
hold on;

%Plot Arrows for first Step
plotRef = gobjects(1,2*n+1);
protRef_circle = gobjects(1,2*n+1);
plotRef_arrowtips = gobjects(1,2*n+1);
plotRef(1) = plot([Vector_position0(1,1,1),Vector_position0(1,3,1)],[Vector_position0(1,2,1),Vector_position0(1,4,1)],'Color',[1,1,1]);
r = calcVectorLength(Vector_position0(1,1,1),Vector_position0(1,2,1),Vector_position0(1,3,1),Vector_position0(1,4,1));
protRef_circle(1) = circle(Vector_position0(1,1,1),Vector_position0(1,2,1),r);
plotRef_arrowtips(1) = drawArrowTip(Vector_position0(1,:,1));
for i=1:n
    plotRef(2*i) = plot([Vector_position_positive(i,1,1),Vector_position_positive(i,3,1)],[Vector_position_positive(i,2,1),Vector_position_positive(i,4,1)],'Color',[1,1,1]);
    r = calcVectorLength(Vector_position_positive(i,1,1),Vector_position_positive(i,2,1),Vector_position_positive(i,3,1),Vector_position_positive(i,4,1));
    protRef_circle(2*i) = circle(Vector_position_positive(i,1,1),Vector_position_positive(i,2,1),r);
    plotRef_arrowtips(2*i) = drawArrowTip(Vector_position_positive(i,:,1));

    plotRef(2*i+1) = plot([Vector_position_negative(i,1,1),Vector_position_negative(i,3,1)],[Vector_position_negative(i,2,1),Vector_position_negative(i,4,1)],'Color',[1,1,1]);
    r = calcVectorLength(Vector_position_negative(i,1,1),Vector_position_negative(i,2,1),Vector_position_negative(i,3,1),Vector_position_negative(i,4,1));
    protRef_circle(2*i+1) = circle(Vector_position_negative(i,1,1),Vector_position_negative(i,2,1),r);
    plotRef_arrowtips(2*i+1) = drawArrowTip(Vector_position_negative(i,:,1));
end
hold off;

%Update plot for animation
for aa = 1:k
    plotRef(1).XData = [Vector_position0(1,1,aa),Vector_position0(1,3,aa)];
    plotRef(1).YData = [Vector_position0(1,2,aa),Vector_position0(1,4,aa)];
    updateArrowTip(Vector_position0(1,:,aa),plotRef_arrowtips(1),7);
    pointer = 2;
    for i=1:n
        plotRef(pointer).XData = [Vector_position_positive(i,1,aa),Vector_position_positive(i,3,aa)];
        plotRef(pointer).YData = [Vector_position_positive(i,2,aa),Vector_position_positive(i,4,aa)];
        r = calcVectorLength(Vector_position_positive(i,1,aa),Vector_position_positive(i,2,aa),Vector_position_positive(i,3,aa),Vector_position_positive(i,4,aa));
        updateCirclePos(protRef_circle(pointer),Vector_position_positive(i,1,aa),Vector_position_positive(i,2,aa),r)
        updateArrowTip(Vector_position_positive(i,:,aa),plotRef_arrowtips(2*i),7/sqrt(i));
        pointer = pointer + 1;
        
        plotRef(pointer).XData = [Vector_position_negative(i,1,aa),Vector_position_negative(i,3,aa)];
        plotRef(pointer).YData = [Vector_position_negative(i,2,aa),Vector_position_negative(i,4,aa)];
        r = calcVectorLength(Vector_position_negative(i,1,aa),Vector_position_negative(i,2,aa),Vector_position_negative(i,3,aa),Vector_position_negative(i,4,aa));
        updateCirclePos(protRef_circle(pointer),Vector_position_negative(i,1,aa),Vector_position_negative(i,2,aa),r)
        updateArrowTip(Vector_position_negative(i,:,aa),plotRef_arrowtips(2*i+1),7/sqrt(i));
        pointer = pointer + 1;
        
    end
    functionH.XData = f(1:aa,1);
    functionH.YData = f(1:aa,2);
    %xlim([Vector_position_negative(n,3,aa)-25 Vector_position_negative(n,3,aa)+25])
    %ylim([Vector_position_negative(n,4,aa)-25 Vector_position_negative(n,4,aa)+25])
    xlim([-150 150]);ylim([-150 150]);
    drawnow;
end

function [Cn_real, Cn_imaginary] = calc_Cn(data_points, j)
    noOfPoints = size(data_points,1);
    dt = 1/(noOfPoints-1);
    t = 0;
    Cn_real = (data_points(1,1)*cos(-j*2*pi*t) - data_points(1,2)*sin(-j*2*pi*t))*dt/2;
    Cn_imaginary = (data_points(1,2)*cos(-j*2*pi*t) + data_points(1,1)*sin(-j*2*pi*t))*dt/2;
    
    for i=2:noOfPoints-1
        t = t + dt;
        Cn_real = Cn_real + (data_points(i,1)*cos(-j*2*pi*t) - data_points(i,2)*sin(-j*2*pi*t))*dt;
        Cn_imaginary = Cn_imaginary + (data_points(i,2)*cos(-j*2*pi*t) + data_points(i,1)*sin(-j*2*pi*t))*dt;
    end
    t = t + dt;
    Cn_real = Cn_real + (data_points(end,1)*cos(-j*2*pi*t) - data_points(end,2)*sin(-j*2*pi*t))*dt/2;
    Cn_imaginary = Cn_imaginary + (data_points(end,2)*cos(-j*2*pi*t) + data_points(end,1)*sin(-j*2*pi*t))*dt/2;

end

function term = getFourierRealTerm(C, j, t)
    term = C(1)*cos(j*2*pi*t) - C(2)*sin(j*2*pi*t);
end


function term = getFourierImagTerm(C, j, t)
    term = C(2)*cos(j*2*pi*t) + C(1)*sin(j*2*pi*t);
end

function h = circle(x,y,r)
    d = r*2;
    px = x-r;
    py = y-r;
    h = rectangle('Position',[px py d d],'Curvature',[1,1],"EdgeColor",[0.7,0.7,0.7],"LineWidth",0.3,"LineStyle",":");
end

function updateCirclePos(h, x, y, r)
    d = r*2;
    px = x-r;
    py = y-r;
    h.Position = [px py d d];
end

function l = calcVectorLength(x1,y1,x2,y2)
    l = sqrt((x1-x2)^2+(y1-y2)^2);
end

function h = drawArrowTip(vector)
    theta = atan((vector(4) - vector(2))/ (vector(3) - vector(1)));
    l = 7;
    h = polyshape([vector(3) vector(3)-l*cos(theta-pi/6) vector(3)-l*cos(theta+pi/6)],[vector(4) vector(4)-l*sin(theta-pi/6) vector(4)-l*sin(theta+pi/6)]);
    h = plot(h, "FaceColor","white",'FaceAlpha',1);
end

function updateArrowTip(vector,h,l)
    l = calcVectorLength(vector(1),vector(2),vector(3),vector(4));
    l = l/3;
    l = min(l,7);
    theta = atan((vector(4) - vector(2))/ (vector(3) - vector(1)));
    if (vector(4) - vector(2)) > 0 && (vector(3) - vector(1)) > 0
        quad = 1;
    elseif (vector(4) - vector(2)) > 0 && (vector(3) - vector(1)) < 0
        quad = 2;
    elseif (vector(4) - vector(2)) < 0 && (vector(3) - vector(1)) < 0
        quad = 3;
    elseif (vector(4) - vector(2)) < 0 && (vector(3) - vector(1)) > 0
        quad = 4;
    else
        quad = 1;
    end
    if quad == 1
        pgon = polyshape([vector(3) vector(3)-l*cos(theta-pi/6) vector(3)-l*cos(theta+pi/6)],[vector(4) vector(4)-l*sin(theta-pi/6) vector(4)-l*sin(theta+pi/6)]);
    elseif quad == 2
        pgon = polyshape([vector(3) vector(3)+l*cos(theta-pi/6) vector(3)+l*cos(theta+pi/6)],[vector(4) vector(4)+l*sin(theta-pi/6) vector(4)+l*sin(theta+pi/6)]);
    elseif quad == 3
        pgon = polyshape([vector(3) vector(3)+l*cos(theta-pi/6) vector(3)+l*cos(theta+pi/6)],[vector(4) vector(4)+l*sin(theta-pi/6) vector(4)+l*sin(theta+pi/6)]);
    else
        pgon = polyshape([vector(3) vector(3)-l*cos(theta-pi/6) vector(3)-l*cos(theta+pi/6)],[vector(4) vector(4)-l*sin(theta-pi/6) vector(4)-l*sin(theta+pi/6)]);
    end
    h.Shape = pgon;
    %h.Shape.Vertices = [vector(3) vector(3)-l*cos(theta-pi/6) vector(3)-l*cos(theta+pi/6);vector(4) vector(4)-l*sin(theta-pi/6) vector(4)-l*sin(theta+pi/6)]';
end