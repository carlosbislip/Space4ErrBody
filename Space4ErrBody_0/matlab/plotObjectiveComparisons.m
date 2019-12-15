function [  ] = plotObjectiveComparisons( compilation )
disp('   Objective Comparisons')

mkdir(strcat(compilation(1).workingFolderPath,'figures/'));


%%
for p = 1:numel(compilation)
    
    %%
    extremesAndConstraintsFieldNames  = fieldnames(compilation(p).evolutions(end).population(end).extremesAndConstraints);
    
    X = nan(compilation(p).evolutions(1).population(1).size.collective,length(extremesAndConstraintsFieldNames),numel(compilation(p).evolutions));
    championVector = nan(1,length(extremesAndConstraintsFieldNames),numel(compilation(p).evolutions));
    championMagnitude = nan(1,length(extremesAndConstraintsFieldNames),numel(compilation(p).evolutions));
    fillerZ = ones(compilation(p).evolutions(1).population(1).size.collective,1);
    maxX = nan( numel(compilation(p).evolutions),length(extremesAndConstraintsFieldNames));
    minX = nan( numel(compilation(p).evolutions),length(extremesAndConstraintsFieldNames));
    
    
    for k = 1:numel(compilation(p).evolutions)
        for ii = 1:compilation(p).evolutions(k).population(1).size.collective
            if compilation(p).evolutions(k).population(ii).rankVector == 1
                for j = 1:length(extremesAndConstraintsFieldNames)
                    championVector(1,j,k) = compilation(p).evolutions(k).population(ii).extremesAndConstraints.(extremesAndConstraintsFieldNames{j}).value;
                end
            end
            if compilation(p).evolutions(k).population(ii).rankMagnitude == 1
                
                for j = 1:length(extremesAndConstraintsFieldNames)
                    championMagnitude(1,j,k) = compilation(p).evolutions(k).population(ii).extremesAndConstraints.(extremesAndConstraintsFieldNames{j}).value;
                end
            end
            for j = 1:length(extremesAndConstraintsFieldNames)
                X(ii,j,k) = compilation(p).evolutions(k).population(ii).extremesAndConstraints.(extremesAndConstraintsFieldNames{j}).value;
            end
        end
    end
    
    for k = 1:numel(compilation(p).evolutions)
        maxX(k,:) = max(X(:,:,k));
        minX(k,:) = min(X(:,:,k));
    end
    
    championVector = reshape(permute(championVector,[1 3 2]),[],size(championVector,2));
    championMagnitude = reshape(permute(championMagnitude,[1 3 2]),[],size(championMagnitude,2));
    
    
    kk = 0;
    
    
    %%
    if numel(compilation(p).evolutions) > 1
        kk = kk + 1;
        for j = 1:length(extremesAndConstraintsFieldNames)
            kk = kk + 1;
            
            fig_num = p*100 + 3466000 + kk;
            figure(fig_num)
            %set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set(figure(fig_num),'units','pixels','position',[2100,0,560,420])
            set (gca,'Fontsize',15)
            
            field = extremesAndConstraintsFieldNames{j};
            variableLabel = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field).variableLabel;
            figureSaveName = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field).figureSaveNameContent;
            units = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field).units;
            limits = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field).limits;
            scalingFactor = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field).scalingFactor;
            ylim(limits)
            if numel(compilation(p).evolutions) > 1
                max_evolutions = 25*ceil(max((numel(compilation(p).evolutions)-1))/25);
                xlim([0 max_evolutions])
                set(gca,'XTick', round(0:max_evolutions/10:max_evolutions));
            end
            %title({strcat(strrep(convertCharsToStrings(compilation(p).case),'_',' '));strcat('Generational History:_{ }',variableLabel)})
            
            xlabel('Generation','Interpreter','latex') % x-axis label
            ylabel(strcat(variableLabel,{' '},units),'Interpreter','latex') % y-axis label
            %zlabel(extremesAndConstraintsFieldNames{combinations(kk,2)}) % z-axis label
            
            hold on
            grid on
            
            m={'+', 'diamond', 'o', '*', '.', 'x', 'square', 'v', '^', '>', '<', 'pentagram'};
            scatter(linspace(0,length(compilation(p).evolutions),length(compilation(p).evolutions)),maxX(:,j)/scalingFactor,'Marker',m{1})
            scatter(linspace(0,length(compilation(p).evolutions),length(compilation(p).evolutions)),minX(:,j)/scalingFactor,'Marker',m{4})
            scatter(linspace(0,length(compilation(p).evolutions),length(compilation(p).evolutions)),championVector(:,j)/scalingFactor,'Marker',m{2},'MarkerFaceColor', 'flat')
            scatter(linspace(0,length(compilation(p).evolutions),length(compilation(p).evolutions)),championMagnitude(:,j)/scalingFactor,'Marker',m{end},'MarkerFaceColor', 'flat')
            
            hold off
            
            
            legendtext = [{'Maximum'};{'Minimum'};{'Vector Champion'};{'Magnitude Champion'}];
            legend(legendtext,'Location','northwest')
            
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'objectiveGenerationalHistory_',...
                figureSaveName,...
                '_Case',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            close(fig_num);
            
        end
        
    end
    %%
    
    combinations = ...
        [1 8;...
        1 9;...
        1 10;...
        1 11;...
        1 12;...
        1 21;...
        1 22;...
        1 23;...
        1 24];
    
    
    for kk = 1:size(combinations,1)
        
        fig_num = p*100 + 3566000 + kk;
        figure(fig_num)
        %set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
        %set(figure(fig_num),'units','pixels','position',[2100,0,560,420])
        set (gca,'Fontsize',15)
        
        Xx = combinations(kk,1);
        Xy = combinations(kk,2);
        field_x = extremesAndConstraintsFieldNames{Xx};
        variableLabel_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).variableLabel;
        figureSaveName_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).figureSaveNameContent;
        units_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).units;
        limits_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).limits;
        tick_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).tick;
        scalingFactor_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).scalingFactor;
        
        field_y = extremesAndConstraintsFieldNames{Xy};
        variableLabel_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).variableLabel;
        figureSaveName_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).figureSaveNameContent;
        units_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).units;
        limits_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).limits;
        tick_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).tick;
        scalingFactor_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).scalingFactor;
        
        
        ylim(limits_y)
        xlim(limits_x)
        xlabel(string(strcat(variableLabel_x,{' '},units_x)),'Interpreter','latex') % x-axis label
        ylabel(string(strcat(variableLabel_y,{' '},units_y)),'Interpreter','latex') % y-axis label
        set(gca,'TickLabelInterpreter','latex')
        
        set(gca,'YTick', round(limits_y(1):tick_y:limits_y(2)));
        
        %title({strcat(strrep(convertCharsToStrings(compilation(p).case),'_',' '));string(strcat(variableLabel_y,{' vs '},variableLabel_x))})
        
        hold on
        grid on
        
        m = {'+', 'diamond', 'o', '*','square', '<','.', 'x', 'v', '^', '>', 'pentagram'};
        
        %XX = reshape(permute(X,[1 3 2]),[],size(X,2));
        q = 1;
        if numel(compilation(p).evolutions) == 1
            
            for i = 1%[1 size(X,3)]
                g(q) = scatter(X(:,Xx,i)/scalingFactor_x,X(:,Xy,i)/scalingFactor_y,'Marker',m{q});
                q = q + 1;
            end
            for i = [1 size(X,3)]
                % g(q) = scatter(championVector(i,Xx,:)/scalingFactor_x,championVector(i,Xy,:)/scalingFactor_y,'Marker',m{q},'MarkerFaceColor', 'flat');
                %  scatter(championMagnitude(i,Xx,:)/scalingFactor_x,championMagnitude(i,Xy,:)/scalingFactor_y,'Marker',m{q+1},'MarkerFaceColor', 'flat')
                q = q + 2;
            end
        else
            
            for i = [1 size(X,3)]
                % scatter(XX(i,Xx,:)/scalingFactor_x,XX(i,Xy,:)/scalingFactor_y,'Marker',m{1})
                g(q) = scatter(X(:,Xx,i)/scalingFactor_x,X(:,Xy,i)/scalingFactor_y,'Marker',m{q});
                q = q + 1;
            end
            for i = [1 size(X,3)]
                % g(q) = scatter(championVector(i,Xx,:)/scalingFactor_x,championVector(i,Xy,:)/scalingFactor_y,'Marker',m{q},'MarkerFaceColor', 'flat');
                %  scatter(championMagnitude(i,Xx,:)/scalingFactor_x,championMagnitude(i,Xy,:)/scalingFactor_y,'Marker',m{q+1},'MarkerFaceColor', 'flat')
                q = q + 2;
            end
        end
        
        hold off
        
        axP = get(gca,'Position');
        set(gca, 'Position', axP)
        %legendText = [
        %   strcat({'Generation: 0'});
        %   strcat({'Generation: '},num2str(numel(compilation(p).evolutions)));
        %   {'Vector Champion - Generation: 0 '};
        %   {'Magnitude Champion - Generation: 0 '};
        %   strcat({'Vector Champion - '},{'Generation: '},num2str(numel(compilation(p).evolutions)));
        %   strcat({'Magnitude Champion - '},{'Generation: '},num2str(numel(compilation(p).evolutions)))];
        %legend(legendtext,'Location','northwest')
        %[~,object_h,~,~] = legendflex(g,legendText,'Interpreter','latex');
        
        saveas(...
            figure(fig_num),...
            strcat(...
            compilation(p).figurePath,...
            'objectiveComparisons_',...
            strcat(figureSaveName_x,'_vs_',figureSaveName_y),...
            '_Case',...
            convertCharsToStrings(compilation(p).case),...
            '.png'),...
            'png');
        close(fig_num);
        
    end
    
    %%
    
    combinations = ...
        [1 9 8;...
        1 10 8;...
        1 11 8;...
        1 12 8;...
        1 21 8;...
        1 22 8;...
        1 23 8;...
        1 24 8 ];
    
    if compilation(1).evolutions(1).population(1).size.collective < 10000
        
        
        for kk = 1:size(combinations,1)
            
            fig_num = p*100 + 3566000 + kk;
            figure(fig_num)
            %set(figure(fig_num),'units','pixels','position',[0,0,1200,600])
            set(figure(fig_num),'units','pixels','position',[2100,0,560,420])
            set (gca,'Fontsize',15)
            
            Xx = combinations(kk,1);
            Xy = combinations(kk,2);
            Xz = combinations(kk,3);
            field_x = extremesAndConstraintsFieldNames{Xx};
            variableLabel_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).variableLabel;
            figureSaveName_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).figureSaveNameContent;
            units_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).units;
            limits_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).limits;
            scalingFactor_x = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_x).scalingFactor;
            
            field_y = extremesAndConstraintsFieldNames{Xy};
            variableLabel_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).variableLabel;
            figureSaveName_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).figureSaveNameContent;
            units_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).units;
            limits_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).limits;
            scalingFactor_y = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_y).scalingFactor;
            
            field_z = extremesAndConstraintsFieldNames{Xz};
            variableLabel_z = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_z).variableLabel;
            figureSaveName_z = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_z).figureSaveNameContent;
            units_z = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_z).units;
            limits_z = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_z).limits;
            scalingFactor_z = compilation(p).evolutions(end).population(end).extremesAndConstraints.(field_z).scalingFactor;
            
            
            ylim(limits_y)
            xlim(limits_x)
            zlim(limits_z)
            
            xlabel(string(strcat(variableLabel_x,{' '},units_x)),'Interpreter','latex') % x-axis label
            ylabel(string(strcat(variableLabel_y,{' '},units_y)),'Interpreter','latex') % y-axis label
            zlabel(string(strcat(variableLabel_z,{' '},units_z)),'Interpreter','latex') % z-axis label
            set(gca,'TickLabelInterpreter','latex')
            
            
            %title({string(strcat(strrep(convertCharsToStrings(compilation(p).case),'_',' ')));string(strcat(variableLabel_y,{' vs '},variableLabel_x,{' vs '},variableLabel_z))})
            
            hold on
            grid on
            
            m = {'+', 'diamond', 'o', '*','square', '<','.', 'x', 'v', '^', '>', 'pentagram'};
            
            i = 1;
            for k = [1 size(X,3)]
                scatter3(X(:,Xx,k)/scalingFactor_x,X(:,Xy,k)/scalingFactor_y,X(:,Xz,k)/scalingFactor_z,'Marker',m{i})
                i = i + 1;
            end
            
            i = 1;
            for k = [1 size(X,3)]
                
                scatter3(X(:,Xx,k)/scalingFactor_x,limits_y(2)*fillerZ,X(:,Xz,k)/scalingFactor_z,'MarkerFacecolor',[0.7 0.7 0.7],'MarkerEdgecolor',[0.7 0.7 0.7],'Marker',m{i})
                scatter3(limits_x(1)*fillerZ,X(:,Xy,k)/scalingFactor_y,X(:,Xz,k)/scalingFactor_z,'MarkerFacecolor',[0.7 0.7 0.7],'MarkerEdgecolor',[0.7 0.7 0.7],'Marker',m{i})
                scatter(X(:,Xx,k)/scalingFactor_x,X(:,Xy,k)/scalingFactor_y,'MarkerFacecolor',[0.7 0.7 0.7],'MarkerEdgecolor',[0.7 0.7 0.7],'Marker',m{i})
                
                i = i + 1;
            end
            
            plot3([0 0],[limits_y(2) limits_y(2)],[limits_z(1) limits_z(2)],'Color','k', 'LineStyle','-')
            plot3([0 0],[limits_y(1) limits_y(2)],[0 0],'Color','k', 'LineStyle','-')
            plot3([limits_x(1) limits_x(2)],[limits_y(2) limits_y(2)],[0 0],'Color','k', 'LineStyle','-')
            
            
            hold off
            
            axP = get(gca,'Position');
            set(gca, 'Position', axP)
            % legendtext = [
            %     strcat({'Generation: 0'});
            %     strcat({'Generation: '},num2str(numel(compilation(p).evolutions)))];
            %legendText = [{'Ascent'};{'Violation during ascent'};{'Descent'};{'Violation during descent'}];
            %[~,object_h,~,~] = legendflex(g,legendText,'Interpreter','latex');
            
            %legend(legendtext,'Location','northwest')
            %view([13 49])
            view([50 17])
            
            saveas(...
                figure(fig_num),...
                strcat(...
                compilation(p).figurePath,...
                'objectiveComparisons_3D_',...
                strcat(figureSaveName_x,'_vs_',figureSaveName_z,'_vs_',figureSaveName_y),...
                '_Case',...
                convertCharsToStrings(compilation(p).case),...
                '.png'),...
                'png');
            close(fig_num);
        end
        
        
    end
    
    
    
    
end