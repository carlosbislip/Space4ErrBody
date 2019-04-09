%% Trim Analysis of


[X,Y] = meshgrid(Mach_b,alpha_b);





%  for k = numel(compilation(p).evolutions):numel(compilation(p).evolutions)
fig_num = 3452345;
figure(fig_num)
set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
set (gca,'Fontsize',15)
% title(strcat('Mach through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
ylim([0 pi/4])
xlim([0 20])
xlabel('Propagation Time (s)') % x-axis label
ylabel('Mach (-)') % y-axis label
set(gca,'YTick', 0:pi/16:pi/4);
set(gca,'XTick', 0:5:20);
hold on

surf(X,Y,Cm0)
view(3)
grid on
hold off
%     saveas(...
%         figure(fig_num),...
%         strcat(...
%         compilation(p).mainpath,...
%         '/figures/mach_v_T_Evolution_',...
%         num2str(k - 1),...
%         '_Set',...
%         convertCharsToStrings(compilation(p).set),...
%         '.png'),...
%         'png');
%    close(fig_num);

[Xb,Yb] = meshgrid(Mach_cs,alpha_cs,db_cs);
[Xb_grid,Yb_grid,Zb_grid] = ndgrid(Mach_cs,alpha_cs',db_cs);
XI = [Xb_grid(:) Yb_grid(:) Zb_grid(:)];

V = dCmb+Cm0;
V = dCmb;
ii = 1;
jj = 1;
kk = 1;
close all
clear Cm_int_alpha Cm_int_bodyflap Cm_int_total
alpha    = [0 5 10 15 20 25 30 35 40 45];
mach     = [1.2 1.5 2 3 5 10 20];
bodyflap = [-20 -10 0 10 20 30];
elevon   = [-40 -30 -20 -10 0 10 20 30 40];
Cm_int_0 = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
Cm_int_bodyflap = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
Cm_int_elevon = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));
Cm_int_matrix = zeros(numel(alpha),numel(mach),numel(bodyflap),numel(elevon));

for i = 1:numel(alpha)
    
    xq1 = deg2rad(alpha(i));
    
    for j = 1:numel(mach)
        
        xq2 = mach(j);
        Cm0_int = interpn(alpha_cs',Mach_cs,Cm0,xq1,xq2,'spline');
        
        for k = 1:numel(bodyflap)
            
            xq31 = deg2rad(bodyflap(k));
            dCmb_int = interpn(alpha_cs',Mach_cs,db_cs,dCmb,xq1,xq2,xq31,'spline');
            
            for l = 1:numel(elevon)
                
                xq32 = deg2rad(elevon(l));
                dCme_int = interpn(alpha_cs',Mach_cs,dw_cs,dCmw,xq1,xq2,xq32,'spline');
                
                Cm_int_0(i,j,k,l) = Cm0_int;
                Cm_int_bodyflap(i,j,k,l) = dCmb_int;
                Cm_int_elevon(i,j,k,l) = dCme_int;
                Cm_int_matrix(i,j,k,l) = Cm0_int + dCmb_int + dCme_int;
                
            end
        end
    end
end

machrange = linspace(1.2,20,(20-1.2)/.25);
alpharange = deg2rad(linspace(0,45,(45-0)/.25));
bodyflaprange = deg2rad(linspace(-20,30,(30--20)/.25));
elevonrange = deg2rad(linspace(-40,40,(40--40)/.5));

pair1 = [ ];
pair2 = [ ];
pair3 = [ ];
p = 0;
for i = 1:numel(machrange)
    
    %
    %     if norm(pair1(1,:)) > 0
    %         p = size(pair1(p,:),1);
    %     else
    %         p = p + 1;
    %     end
    %
    for j = 1:numel(alpharange)
        
        Cm_int_total = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,machrange(i),alpharange(j),0.0,0.0,'spline');
        %  machrange(i)
        %  rad2deg(alpharange(j))
        if abs(Cm_int_total) < 0.0001
            pair1 = [ pair1 ; machrange(i) rad2deg(alpharange(j)) Cm_int_total ];
            %             p = p + 1;
        else
            %             if norm(pair2(1,:)) > 0
            %                 p = size(pair2(p,:),1);
            %             else
            %              %   p = p + 1;
            %             end
            for k = 1:numel(bodyflaprange)
                
                % rad2deg(bodyflaprange(k))
                
                Cm_int_total = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,machrange(i),alpharange(j),bodyflaprange(k),0.0,'spline');
                if abs(Cm_int_total) < 0.0001
                    pair2 = [ pair2; machrange(i) rad2deg(alpharange(j)) rad2deg(bodyflaprange(k)) Cm_int_total ];
                    %   p = p + 1;
                else
                    %                     if norm(pair3(1,:)) > 0
                    %                         p = size(pair3(p,:),1);
                    %                     else
                    %                         p = p + 1;
                    %                     end
                    for l = 1:numel(elevonrange)
                        
                        %rad2deg( elevonrange(l))
                        
                        Cm_int_total = interpn(alpha_cs',Mach_cs,db_cs,dw_cs,Cm_int_matrix,machrange(i),alpharange(j),bodyflaprange(k),elevonrange(l),'spline');
                        if abs(Cm_int_total) < 0.0001
                            pair3 = [pair3 ; machrange(i) rad2deg(alpharange(j)) rad2deg(bodyflaprange(k)) rad2deg(elevonrange(l)) Cm_int_total ]
                        end
                        
                        
                    end
                end
                
            end
        end
    end
end

b = repmat( reshape(a, size(a,1), size(a,2), 1, size(a,3)), 1, 1, size(a,2), 1);


x = -5:0.8:5;
y = x';
z = sin(x.^2 + y.^2) ./ (x.^2 + y.^2);
F = griddedInterpolant(ndgrid({Mach_cs,alpha_cs,db_cs}),V);
xq = -5:0.1:5;
yq = xq';
vq = F({xq,yq});

[X1,X2,X3] = ndgrid({Mach_cs,alpha_cs,db_cs})

fig_num = fig_num + 1;
figure(fig_num)
set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
set (gca,'Fontsize',15)
% title(strcat('Mach through Time - Evolution:_{ }',num2str(k - 1),' - ',strrep(convertCharsToStrings(compilation(p).set),'_',' ')))
ylim([0 pi/4])
xlim([0 20])
xlabel('Propagation Time (s)') % x-axis label
ylabel('Mach (-)') % y-axis label
set(gca,'YTick', 0:pi/16:pi/4);
set(gca,'XTick', 0:5:20);
hold on
%plot3(Xb(:,:,1),Yb(:,:,1),Zb(:,:,2))
for k = 1:6
    surf(Xb(:,:,k),Yb(:,:,k),dCmb(:,:,k))
end
view(3)
grid on
hold off

