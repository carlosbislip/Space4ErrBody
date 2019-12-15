function [  ] = plotAlphaMachForDiscussion( compilation )


%% Angle of Attack vs. Mach Number - per Generation
p = 1;
k = 1;


fig_num = p*100 + k*1 + 100000;
figure(fig_num)
set (gca,'Fontsize',15)

field_x = 'machNumber';
field_x_struct = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_x);
variableLabel_x = field_x_struct.variableLabel;
figureSaveName_x = field_x_struct.figureSaveNameContent;
units_x = field_x_struct.units;
limits_x = field_x_struct.limits;
scalingFactor_x = field_x_struct.scalingFactor;
tick_x = field_x_struct.tick;

field_y = 'angleOfAttack';
field_y_struct = compilation(p).evolutions(end).population(end).dependentVariableTimeHistory.(field_y);
variableLabel_y = field_y_struct.variableLabel;
figureSaveName_y = field_y_struct.figureSaveNameContent;
units_y = field_y_struct.units;
limits_y = field_y_struct.limits;
scalingFactor_y = field_y_struct.scalingFactor;
tick_y = field_y_struct.tick;

xlabel(string(strcat(variableLabel_x,{' '},units_x)),'Interpreter','latex') % x-axis label
ylabel(string(strcat(variableLabel_y,{' '},units_y)),'Interpreter','latex') % y-axis label

ylim(limits_y)
xlim(limits_x)

set(gca,'XTick', limits_x(1):tick_x:limits_x(2));
set(gca,'YTick', limits_y(1):tick_y:limits_y(2));
set(gca,'TickLabelInterpreter','latex')

%set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
%set (gca,'Fontsize',15)
%title(strcat('Angle of Attack vs. Mach Number - Generation:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).case),'_',' ')))
%ylim([0 50])
%xlim([0 30])
%xlabel('Mach Number (-)') % x-axis label
%ylabel('Angle of Attack (deg)') % y-axis label
%set(gca,'YTick', 0:10:50);
%set(gca,'XTick', 0:2.5:30);
hold on
grid on

x = compilation(p).evolutions(k).Common.Bounds.AngleOfAttack.machNumber;
y_ub = compilation(p).evolutions(k).Common.Bounds.AngleOfAttack.upperBound;
y_lb = compilation(p).evolutions(k).Common.Bounds.AngleOfAttack.lowerBound;
plot(x,y_ub,'k','LineWidth',2);
plot(x,y_lb,'k','LineWidth',2);

shadedColor = [ 0.5843 0.8157 0.9882 ];

%[h(numel(h) + 1),msg]=jbfill(transpose(x),transpose(y_ub),transpose(y_lb),shadedColor,shadedColor/2,0,.5);
[h,msg]=jbfill(transpose(x),transpose(y_ub),transpose(y_lb),shadedColor,shadedColor/2,0,.5);

%patches = find( isgraphics(h,'patch') );
hatchAngle = [135 45];
hatchType = [{'single'} {'single'}];

hatchfill2(h,hatchType{1},'HatchAngle',hatchAngle(1),'HatchDensity',40,'HatchColor','k','HatchLineWidth',.5);

hold off
saveas(...
    figure(fig_num),...
    strcat(...
    figurePath,...
    'HORUS_alpha_mach_envelope','.png'),'png');
%     close(fig_num);





end