function [  ] = plot_Validation_Trajectories( compilation )

lon_i_deg = -106.7;
lon_i_rad = deg2rad(-22.37);
lon_i_deg = -22.37;
lat_f_deg = 5;
lon_f_deg = -53;
validation = 1;

xunit = compilation.rawData.targetCircleCoordinates(:,3);
yunit = compilation.rawData.targetCircleCoordinates(:,2);

xmax = max(xunit);
xmin = min(xunit);


%% Validation Trajectory on 2D Map

fig_num = 3100;
figure(fig_num)
set(figure(fig_num),'units','pixels','position',[0,0,1200,600])

ax1 = gca;

hold on
xlabel(string(strcat({'$\tau$ $(\mathrm{deg})$'})),'Interpreter','latex') % x-axis label
ylabel(string(strcat({'$\delta$ $(\mathrm{deg})$'})),'Interpreter','latex') % y-axis label
img = imread('img.jpg');
imagesc([-180 180], [-90 90], (flipud(img)));

limit_x = [-120 -50];
x_width = limit_x(2) - limit_x(1);
y_width = x_width/2;
limit_y(1) = -25;
limit_y(2) = limit_y(1) + y_width;
xlim(limit_x)
ylim(limit_y)

fontSize = 20;
set (gca,'Fontsize',fontSize)

set(gca,'YTick', limit_y(1):5:limit_y(2));
set(gca,'XTick', limit_x(1):x_width/10:limit_x(2));

set(gca,'TickLabelInterpreter','latex')

x = compilation(1).evolutions(1).population(1).dependentVariableTimeHistory.longitudeAngle.value;
y = compilation(1).evolutions(1).population(1).dependentVariableTimeHistory.latitudeAngle.value;
plot(ax1, x,y,'g','LineWidth',2)

plot(ax1, lon_f_deg,lat_f_deg,'MarkerSize',20);
plot(ax1, lon_f_deg*[1 1],90*[-1 1],'k','LineWidth',2);
plot(ax1,180*[-1 1],lat_f_deg*[1 1],'k','LineWidth',2);
plot(ax1, xunit, yunit,'k','LineWidth',2);
scatter(lon_f_deg,lat_f_deg,100,'r','x');

ax2 = axes('Position',[.15 .68 .2 .2]);
hold(ax2,'on');
box(ax2,'on');
imagesc(ax2,[-180 180], [-90 90], (flipud(img)));
plot(ax2,lon_f_deg,lat_f_deg,'MarkerSize',20);
plot(ax2,lon_f_deg*[1 1],90*[-1 1],'k','LineWidth',2);
plot(ax2,180*[-1 1],lat_f_deg*[1 1],'k','LineWidth',2);
plot(ax2,xunit, yunit,'k','LineWidth',2);
scatter(ax2,lon_f_deg,lat_f_deg,100,'r','x');

plot(ax2, x,y,'g','LineWidth',2)

grid(ax2, 'on');
xlim(ax2,[(lon_f_deg - 5) (lon_f_deg + 27)])
ylim(ax2,[(lat_f_deg - 5) (lat_f_deg + 11)])
limits_x_zoom = [-55 -51];


x_width = limits_x_zoom(2) - limits_x_zoom(1);
y_width = x_width/2;
limits_y_zoom(1) = 4;
limits_y_zoom(2) = limits_y_zoom(1) + y_width;
xlim(ax2,limits_x_zoom)
ylim(ax2,limits_y_zoom)

set(ax2,'Fontsize',fontSize)
set(ax2,'XTick', limits_x_zoom(1):1:limits_x_zoom(2));
set(ax2,'YTick', limits_y_zoom(1):1:limits_y_zoom(2));
set(ax2,'TickLabelInterpreter','latex')
ax2.YAxisLocation = 'right';

hold off

saveName = string(strcat(...
    compilation.figurePath,...
    'validation_trajectory_2D.png'));

saveas(figure(fig_num),saveName,'png');

close(fig_num);


end