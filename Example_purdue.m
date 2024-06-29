clc; clear; close all;

dt = 0.0005;
n = 20;
data_points = [-265.344633726022,-141.490493544651;
-262.436421165635,-131.796451676696;
-257.589400231658,-120.16360143515;
-251.772975110885,-107.561347006809;
-245.956549990112,-94.9590925784665;
-241.109529056135,-82.3568381501255;
-228.507274627793,-82.3568381501255;
-196.516936463542,-82.3568381501255;
-177.128852727632,-83.3262423369215;
-167.434810859677,-83.3262423369215;
-161.618385738904,-73.6322004689665;
-156.771364804927,-61.9993502274205;
-149.985535497358,-45.5194790518965;
-138.352685255813,-15.4679492612366;
-123.81162245388,21.3694098369924;
-107.331751278357,56.2679605616295;
-101.515326157584,70.8090233635625;
-96.6683052236061,82.4418736051075;
-93.7600926632202,92.1359154730625;
-104.42353871797,93.1053196598585;
-115.086984772721,92.1359154730625;
-130.597451761449,93.1053196598585;
-151.924343870949,92.1359154730625;
-167.434810859677,91.1665112862675;
-163.557194112495,100.860553154222;
-158.710173178518,112.493403395768;
-149.985535497358,129.942678758087;
-144.169110376586,146.42254993361;
-138.352685255813,151.269570867587;
-128.658643387858,150.300166680793;
-109.270559651948,150.300166680793;
-84.0660507952651,152.238975054383;
-14.2689493459902,150.300166680793;
52.6199395428989,150.300166680793;
104.967765629856,151.269570867587;
155.376783343221,151.269570867587;
216.449247111337,151.269570867587;
267.827669011498,148.361358307201;
290.123965307795,145.453145746814;
305.634432296522,139.636720626041;
318.236686724864,130.912082944882;
325.022516032432,117.340424329746;
324.053111845637,99.8911489674265;
319.206090911659,79.5336610447215;
308.542644856909,56.2679605616295;
290.123965307795,17.4917930898104;
276.552306692658,2.95073028787746;
259.103031330339,-11.5903325140546;
240.684351781224,-19.3455660084186;
220.326863858519,-27.1007995027825;
193.183546628245,-30.9784162499645;
171.856654518745,-31.9478204367595;
153.43797496963,-30.9784162499645;
120.478232618584,-31.9478204367595;
97.2125321354919,-30.9784162499645;
61.3445772240589,-30.9784162499645;
28.3848348730119,-30.9784162499645;
5.11913438991985,-31.9478204367595;
-6.51371585162616,-31.9478204367595;
-4.57490747803513,-21.2843743820095;
1.24151764273785,-7.71271576687255;
12.8743678842839,26.2164307709695;
24.5072181258299,27.1858349577655;
36.1400683673759,27.1858349577655;
57.4669604768769,27.1858349577655;
78.7938525863768,26.2164307709695;
106.906574003447,27.1858349577655;
143.743933101675,27.1858349577655;
143.743933101675,27.1858349577655;
163.132016837585,33.0022600785375;
178.642483826313,37.8492810125155;
194.152950815041,55.2985563748345;
198.999971749018,79.5336610447215;
187.367121507473,89.2277029126765;
167.979037771563,92.1359154730625;
131.141678673334,92.1359154730625;
94.3043195751048,92.1359154730625;
55.5281521032858,91.1665112862675;
38.0788767409668,92.1359154730625;
34.2012599937848,82.4418736051075;
29.3542390598068,69.8396191767664;
24.5072181258299,58.2067689352205;
16.7519846314659,40.7574935729014;
5.11913438991985,22.3388140237875;
-0.697290730853126,7.79775122185546;
-6.51371585162616,-9.65152414046352;
-14.2689493459902,-24.1925869423966;
-18.1465660931722,-42.6112664915105;
-24.9323954007402,-58.1217334802385;
-31.7182247083092,-71.6933920953755;
-33.6570330819001,-83.3262423369215;
-16.2077577195811,-83.3262423369215;
-0.697290730853126,-83.3262423369215;
21.5990055654428,-82.3568381501255;
40.0176851145579,-83.3262423369215;
34.2012599937848,-96.8979009520575;
24.5072181258299,-116.285984687968;
18.6907930050568,-132.765855863492;
13.8437720710799,-142.459897731446;
-2.63609910444416,-142.459897731446;
-22.0241828403542,-142.459897731446;
-41.4122665762632,-141.490493544651;
-62.7391586857642,-143.429301918242;
-62.7391586857642,-143.429301918242;
-86.9742633556522,-142.459897731446;
-115.086984772721,-142.459897731446;
-134.475068508631,-142.459897731446;
-153.86315224454,-142.459897731446;
-179.067661101223,-143.429301918242;
-200.394553210724,-142.459897731446;
-221.721445320225,-142.459897731446;
-239.170720682544,-141.490493544651;
-253.711783484476,-141.490493544651;
-265.344633726022,-141.490493544651];

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
xlim([-300 400]);ylim([-300 300]);
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
    xlim([-300 400]);ylim([-300 300]);
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
    h = rectangle('Position',[px py d d],'Curvature',[1,1],"EdgeColor",[0.7,0.7,0.7],"LineWidth",0.1,"LineStyle",":");
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
    l = min(l,20);
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