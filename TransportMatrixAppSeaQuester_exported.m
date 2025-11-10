classdef TransportMatrixAppSeaQuester_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        gCm2yearLabel                 matlab.ui.control.Label
        Image                         matlab.ui.control.Image
        km2Label                      matlab.ui.control.Label
        PgCLabel                      matlab.ui.control.Label
        yearLabel                     matlab.ui.control.Label
        PgCyearLabel                  matlab.ui.control.Label
        secLabel                      matlab.ui.control.Label
        AreaEditField                 matlab.ui.control.NumericEditField
        AreaEditFieldLabel            matlab.ui.control.Label
        p2EditField                   matlab.ui.control.NumericEditField
        p2EditFieldLabel              matlab.ui.control.Label
        p1EditField                   matlab.ui.control.NumericEditField
        p1EditFieldLabel              matlab.ui.control.Label
        VerticalProfileDropDown       matlab.ui.control.DropDown
        VerticalProfileDropDownLabel  matlab.ui.control.Label
        CarbonExportEditField         matlab.ui.control.NumericEditField
        CarbonExportEditFieldLabel    matlab.ui.control.Label
        SequesteredCarbonEditField    matlab.ui.control.NumericEditField
        SequesteredCarbonLabel        matlab.ui.control.Label
        ResidenceTimeEditField        matlab.ui.control.NumericEditField
        ResidenceTimeEditFieldLabel   matlab.ui.control.Label
        RunTimeEditField              matlab.ui.control.NumericEditField
        RunTimeEditFieldLabel         matlab.ui.control.Label
        GoButton                      matlab.ui.control.Button
        dlatEditField                 matlab.ui.control.NumericEditField
        dlatEditFieldLabel            matlab.ui.control.Label
        dlonEditField                 matlab.ui.control.NumericEditField
        dlonEditFieldLabel            matlab.ui.control.Label
        ClearButton                   matlab.ui.control.Button
        FluxEditField                 matlab.ui.control.NumericEditField
        FluxEditFieldLabel            matlab.ui.control.Label
        slatEditField                 matlab.ui.control.NumericEditField
        slatEditFieldLabel            matlab.ui.control.Label
        slonEditField                 matlab.ui.control.NumericEditField
        slonEditFieldLabel            matlab.ui.control.Label
        ProfileAxes                   matlab.ui.control.UIAxes
        SeqMapAxes                    matlab.ui.control.UIAxes
        FluxMapAxes                   matlab.ui.control.UIAxes
        SelectionMapAxes              matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        g % grid data
        grid % grid struct
        TMM % TMM struct
        bathy % bathymetry data
        area % area of grid regions
        dz % vertical depth bin intervals (surface to bottom)
        zstar % depth of the boundaries of depth bins        
        f % flux data matrix (gC / m2 / year)
        f0
        flux = 100 % initial value of flux (gC/m2/year)
        Imax % number of longitude bins
        Jmax % number of latitude bins
        Kmax % number of selected flux sites;
        vmodepara = [0.3,20;1000,100;1000,100;0.1,50]; % initial parameter values for vertical injection profile
        vmodelab = {'b','z0';'z0','zm';'z0','zm';'a','w'}; % initial parameter values for vertical injection profile
        vmodeunits = {'-','m';'m','m';'m','m';'1/s','m/s'}; % initial parameter values for vertical injection profile
        dlat % selection latitude range in indices
        dlon % selection latitude range in indices
    end
    
    methods (Access = public)
        
        function MartinCurve(app)
            
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            load('C:\Users\avis\MATLAB Drive\JeromeSequest\CTL.mat')
            app.TMM = output;
            app.bathy = flipud(sum(app.TMM.M3d,3)); % bottom at app.bathy(bottomindex + 1), bottomindex = 0 => land
            app.area = app.TMM.grid.Areac;
            zw = app.TMM.grid.zw'; % this is the depth of the upper boundaried of depth bins
            app.zstar = [zw; zw(end) + app.TMM.grid.dzt(end)]; % appends the bottom depth
            app.dz = app.zstar(2:end) - app.zstar(1:end-1);
            app.g = app.bathy;
            app.f = flipud(app.TMM.M3d(:,:,1))-1;
            app.f0 = app.f;
            [app.Jmax,app.Imax] = size(app.g);
            imagesc(app.SelectionMapAxes,app.g);
            set(get(app.SelectionMapAxes,'Children'),'HitTest','off'); % Allows mouse clicks on image
            app.SelectionMapAxes.YDir = "reverse"; % doesn't do anything -- ??? -- hence flipud above
            app.SelectionMapAxes.XLim = [1 app.Imax+1]-1/2;
            app.SelectionMapAxes.YLim = [1 app.Jmax+1]-1/2;
            app.FluxEditField.Value = 20;
            imagesc(app.FluxMapAxes,app.f);
            app.FluxMapAxes.XLim = [1 app.Imax+1]-1/2;
            app.FluxMapAxes.YLim = [1 app.Jmax+1]-1/2;
            drawnow;
            app.dlonEditField.Value=2;
            VerticalProfileDropDownValueChanged(app);
    
        end

        % Button down function: SelectionMapAxes
        function SelectionMapAxesButtonDown(app, event)
            point = get(app.SelectionMapAxes,'CurrentPoint');
            isel = round(point(1,1));
            jsel = round(point(1,2));
            app.slonEditField.Value = isel;
            app.slatEditField.Value = jsel;
            app.dlat = app.dlatEditField.Value;
            app.dlon = app.dlonEditField.Value;
            app.flux = app.FluxEditField.Value;
            app.f(jsel,isel)
            if app.f(jsel,isel) > 0
                app.g(jsel,isel) = app.bathy(jsel,isel);
                app.f(jsel,isel) = 0;
            else
                for j = -app.dlat:app.dlat
                    jo = max(1,jsel+j);
                    for i = -app.dlon:app.dlon
                        io = mod(isel+i-1,app.Imax) + 1;
                        if app.bathy(jo,io) > 0
                            app.g(jo,io) = app.flux;
                            app.f(jo,io) = app.flux;
                        end
                    end
                end
            end
            imagesc(app.SelectionMapAxes,app.g);
            imagesc(app.FluxMapAxes,app.f);
            set(get(app.SelectionMapAxes,'Children'),'HitTest','off'); % Allows mouse clicks on image
            drawnow;
            sindex = find(app.f>0)
            app.AreaEditField.Value = sum(app.area(sindex))*1E-6;
        end

        % Button pushed function: ClearButton
        function ClearButtonPushed(app, event)
            app.g = app.bathy;
            app.f = app.f0;
            imagesc(app.FluxMapAxes,app.f);
            imagesc(app.SelectionMapAxes,app.g);
            set(get(app.SelectionMapAxes,'Children'),'HitTest','off'); % Allows mouse clicks on image
            drawnow;
        end

        % Value changed function: dlatEditField
        function dlatEditFieldValueChanged(app, event)
            app.dlat = app.dlatEditField.Value;         
        end

        % Value changed function: dlonEditField
        function dlonEditFieldValueChanged(app, event)
            app.dlon = app.dlonEditField.Value;
        end

        % Button pushed function: GoButton
        function GoButtonPushed(app, event)
            %TR = app.TMM.TR;
            %msk = app.TMM.msk;
            app.f = flipud(app.f); % realign flux matix
            VOL = app.TMM.grid.DXT3d.*app.TMM.grid.DYT3d.*app.TMM.grid.DZT3d;
            V = VOL(app.TMM.msk.pkeep);
            m = size(app.TMM.TR,1);
            sink = zeros(m,1);
            sink(1:length(app.TMM.msk.hkeep)) = 1e10; % instantaneous Surface SINK
            SSINK = spdiags(sink,0,m,m); 
            A = app.TMM.TR-SSINK; % TMM with surface boundary condition
            im = app.VerticalProfileDropDown.ValueIndex;
            % Martin Curve
            z0 = app.p2EditField.Value; b = app.p1EditField.Value;
            app.p1EditFieldLabel.Text = 'b';
            app.p2EditFieldLabel.Text = 'z0';
            %zw = app.TMM.grid.zw'; % this is the depth of the upper boundaried of depth bins
            % zwstar = [zw; zw(end) + app.TMM.grid.dzt(end)]; % appends the bottom depth
            % dz = zwstar(2:end) - zwstar(1:end-1);
            fzstar = (app.zstar/z0).^(-b); % for Martin Curve
            % f = (sin((zo - z)*pi/(2*zm))+1)/2; % for Verical Migration
            % f(z<zo-zm)=1;
            % f(z>zo+zm)=0;
            source = -(fzstar(2:end) - fzstar(1:end-1))./app.dz;
            source(1) = 0; source(end) = fzstar(end-1)./app.dz(end);
            sum(source.*app.dz)
            stairs(app.ProfileAxes,source,-app.zstar(2:end))
            %
            [latindex,lonindex]=find(app.f>0);
            Aimp = zeros(app.Kmax,8)
            app.Kmax = length(lonindex);
            Qz = zeros(size(app.TMM.M3d));
            Qb = zeros(size(app.TMM.M3d));
            AreaSel = 0;
            for k = 1:app.Kmax
                lati = latindex(k);
                loni = lonindex(k);
                ocean = squeeze(app.TMM.M3d(lati,loni,:));
                boti = sum(ocean);
                AreaSel = AreaSel + app.area(lati,loni)
                Aimp(k,:) = [k,lati,loni,boti,app.f(lati,loni),im,app.vmodepara(im,1),app.vmodepara(im,2)];
                if loni < 180
                    %[k, lati,loni,boti]
                    if boti > 0
                        Qz(lati,loni,:) = app.f(lati,loni)*source.*ocean;
                        Qz(lati,loni,boti) = 0;
                        Qb(lati,loni,boti) = app.f(lati,loni)*fzstar(boti)/app.dz(boti);
                        fb(lati,loni) = fzstar(boti);
                    end
                end
            end
            Q = Qz + Qb;
            % i = 11; q = squeeze(Q(latindex(i),lonindex(i),:));
            % figure(7); clf; stairs((q),-zwstar(2:end))
            q_ocim = Q(app.TMM.msk.pkeep);
            q_ocim(isnan(q_ocim)) = 0;
            export = V'*q_ocim / 1e15 % [PgC / yr]
            tic
            cseq = -A \q_ocim;  % Inverts TMM for steady state.
            app.RunTimeEditField.Value = toc;
            TotCseq = V'*cseq / 1e15 % grams to Pg
            SeqTime = TotCseq / export
            app.CarbonExportEditField.Value = export;
            app.SequesteredCarbonEditField.Value = TotCseq;
            app.ResidenceTimeEditField.Value = SeqTime;
            cout = 0*app.TMM.M3d+NaN; % make a 3-d array of NaN
            cout(app.TMM.msk.pkeep) = cseq;
            cout(isnan(cout)) = 0;
            cint = sum(cout.*app.TMM.grid.DZT3d,3); %cc = [cc, cc(:,end)]; % gC/m2
            imagesc(app.SeqMapAxes,flipud(cint));
            app.SeqMapAxes.XLim = [1 app.Imax+1]-1/2;
            app.SeqMapAxes.YLim = [1 app.Jmax+1]-1/2;
            save('TMMout','Aimp')
        end

        % Value changed function: FluxEditField
        function FluxEditFieldValueChanged(app, event)
            app.flux = app.FluxEditField.Value;
            
        end

        % Drop down opening function: VerticalProfileDropDown
        function VerticalProfileDropDownOpening(app, event)
            % vmode = app.VerticalProfileDropDown.ValueIndex
            % app.p1EditField.Value = app.vmodepara(vmode,1);
            % app.p2EditField.Value = app.vmodepara(vmode,2);
            % app.p1EditFieldLabel.Text = app.vmodelab(vmode,1);
            % app.p2EditFieldLabel.Text = app.vmodelab(vmode,2);
        end

        % Value changed function: p1EditField
        function p1EditFieldValueChanged(app, event)
            vmode = app.VerticalProfileDropDown.ValueIndex;
            app.vmodepara(vmode,1) = app.p1EditField.Value;
            
        end

        % Value changed function: p2EditField
        function p2EditFieldValueChanged(app, event)
            vmode = app.VerticalProfileDropDown.ValueIndex;
            app.vmodepara(vmode,2) = app.p2EditField.Value;
            
        end

        % Clicked callback: VerticalProfileDropDown
        function VerticalProfileDropDownClicked(app, event)
            % vmode = app.VerticalProfileDropDown.ValueIndex
            % app.p1EditField.Value = app.vmodepara(vmode,1);
            % app.p2EditField.Value = app.vmodepara(vmode,2);
            % app.p1EditFieldLabel.Text = app.vmodelab(vmode,1);
            % app.p2EditFieldLabel.Text = app.vmodelab(vmode,2);
        end

        % Value changed function: VerticalProfileDropDown
        function VerticalProfileDropDownValueChanged(app, event)
            value = app.VerticalProfileDropDown.Value;
            vmode = app.VerticalProfileDropDown.ValueIndex
            app.p1EditField.Value = app.vmodepara(vmode,1);
            app.p2EditField.Value = app.vmodepara(vmode,2);
            app.p1EditFieldLabel.Text = app.vmodelab(vmode,1);
            app.p2EditFieldLabel.Text = app.vmodelab(vmode,2);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1154 774];
            app.UIFigure.Name = 'MATLAB App';

            % Create SelectionMapAxes
            app.SelectionMapAxes = uiaxes(app.UIFigure);
            title(app.SelectionMapAxes, 'Selection Map')
            app.SelectionMapAxes.YDir = 'reverse';
            app.SelectionMapAxes.ButtonDownFcn = createCallbackFcn(app, @SelectionMapAxesButtonDown, true);
            app.SelectionMapAxes.PickableParts = 'all';
            app.SelectionMapAxes.Position = [14 482 407 273];

            % Create FluxMapAxes
            app.FluxMapAxes = uiaxes(app.UIFigure);
            title(app.FluxMapAxes, 'Flux Map')
            app.FluxMapAxes.Position = [616 482 407 273];

            % Create SeqMapAxes
            app.SeqMapAxes = uiaxes(app.UIFigure);
            title(app.SeqMapAxes, 'Sequestration Map')
            app.SeqMapAxes.Position = [617 23 407 273];

            % Create ProfileAxes
            app.ProfileAxes = uiaxes(app.UIFigure);
            title(app.ProfileAxes, 'Injection  Profile')
            ylabel(app.ProfileAxes, 'depth')
            zlabel(app.ProfileAxes, 'Z')
            app.ProfileAxes.Position = [379 23 215 327];

            % Create slonEditFieldLabel
            app.slonEditFieldLabel = uilabel(app.UIFigure);
            app.slonEditFieldLabel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.slonEditFieldLabel.HorizontalAlignment = 'right';
            app.slonEditFieldLabel.FontWeight = 'bold';
            app.slonEditFieldLabel.FontColor = [0.4667 0.6745 0.1882];
            app.slonEditFieldLabel.Position = [437 707 30 22];
            app.slonEditFieldLabel.Text = 'slon';

            % Create slonEditField
            app.slonEditField = uieditfield(app.UIFigure, 'numeric');
            app.slonEditField.Editable = 'off';
            app.slonEditField.FontWeight = 'bold';
            app.slonEditField.FontColor = [0.4667 0.6745 0.1882];
            app.slonEditField.BackgroundColor = [0.9412 0.9412 0.9412];
            app.slonEditField.Tooltip = {'selected latitude index'};
            app.slonEditField.Position = [482 703 83 22];

            % Create slatEditFieldLabel
            app.slatEditFieldLabel = uilabel(app.UIFigure);
            app.slatEditFieldLabel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.slatEditFieldLabel.HorizontalAlignment = 'right';
            app.slatEditFieldLabel.FontWeight = 'bold';
            app.slatEditFieldLabel.FontColor = [0.4667 0.6745 0.1882];
            app.slatEditFieldLabel.Position = [441 678 26 22];
            app.slatEditFieldLabel.Text = 'slat';

            % Create slatEditField
            app.slatEditField = uieditfield(app.UIFigure, 'numeric');
            app.slatEditField.Editable = 'off';
            app.slatEditField.FontWeight = 'bold';
            app.slatEditField.FontColor = [0.4667 0.6745 0.1882];
            app.slatEditField.BackgroundColor = [0.9412 0.9412 0.9412];
            app.slatEditField.Tooltip = {'selected latitude index'};
            app.slatEditField.Position = [482 676 83 22];

            % Create FluxEditFieldLabel
            app.FluxEditFieldLabel = uilabel(app.UIFigure);
            app.FluxEditFieldLabel.HorizontalAlignment = 'right';
            app.FluxEditFieldLabel.Position = [52 462 28 22];
            app.FluxEditFieldLabel.Text = 'Flux';

            % Create FluxEditField
            app.FluxEditField = uieditfield(app.UIFigure, 'numeric');
            app.FluxEditField.ValueChangedFcn = createCallbackFcn(app, @FluxEditFieldValueChanged, true);
            app.FluxEditField.Tooltip = {'Flux [g C / m^2 / year)'};
            app.FluxEditField.Position = [50 441 69 22];

            % Create ClearButton
            app.ClearButton = uibutton(app.UIFigure, 'push');
            app.ClearButton.ButtonPushedFcn = createCallbackFcn(app, @ClearButtonPushed, true);
            app.ClearButton.BackgroundColor = [0.9412 0.9412 0.9412];
            app.ClearButton.FontSize = 18;
            app.ClearButton.FontWeight = 'bold';
            app.ClearButton.FontColor = [0.851 0.3255 0.098];
            app.ClearButton.Tooltip = {'Refresh stected areas and fluxes'};
            app.ClearButton.Position = [483 623 82 30];
            app.ClearButton.Text = 'Clear';

            % Create dlonEditFieldLabel
            app.dlonEditFieldLabel = uilabel(app.UIFigure);
            app.dlonEditFieldLabel.HorizontalAlignment = 'right';
            app.dlonEditFieldLabel.Position = [224 442 28 22];
            app.dlonEditFieldLabel.Text = 'dlon';

            % Create dlonEditField
            app.dlonEditField = uieditfield(app.UIFigure, 'numeric');
            app.dlonEditField.ValueChangedFcn = createCallbackFcn(app, @dlonEditFieldValueChanged, true);
            app.dlonEditField.Tooltip = {'Selection longitude width (in indices)'};
            app.dlonEditField.Position = [257 441 47 22];
            app.dlonEditField.Value = 1;

            % Create dlatEditFieldLabel
            app.dlatEditFieldLabel = uilabel(app.UIFigure);
            app.dlatEditFieldLabel.HorizontalAlignment = 'right';
            app.dlatEditFieldLabel.Position = [312 442 25 22];
            app.dlatEditFieldLabel.Text = 'dlat';

            % Create dlatEditField
            app.dlatEditField = uieditfield(app.UIFigure, 'numeric');
            app.dlatEditField.ValueChangedFcn = createCallbackFcn(app, @dlatEditFieldValueChanged, true);
            app.dlatEditField.Tooltip = {'Selection latitude width (in indices)'};
            app.dlatEditField.Position = [348 441 47 22];
            app.dlatEditField.Value = 1;

            % Create GoButton
            app.GoButton = uibutton(app.UIFigure, 'push');
            app.GoButton.ButtonPushedFcn = createCallbackFcn(app, @GoButtonPushed, true);
            app.GoButton.FontSize = 18;
            app.GoButton.FontWeight = 'bold';
            app.GoButton.FontColor = [0.851 0.3255 0.098];
            app.GoButton.Tooltip = {'Calculate carbon setquestration '};
            app.GoButton.Position = [484 566 82 30];
            app.GoButton.Text = 'Go';

            % Create RunTimeEditFieldLabel
            app.RunTimeEditFieldLabel = uilabel(app.UIFigure);
            app.RunTimeEditFieldLabel.HorizontalAlignment = 'right';
            app.RunTimeEditFieldLabel.FontWeight = 'bold';
            app.RunTimeEditFieldLabel.Position = [841 432 59 22];
            app.RunTimeEditFieldLabel.Text = 'Run Time';

            % Create RunTimeEditField
            app.RunTimeEditField = uieditfield(app.UIFigure, 'numeric');
            app.RunTimeEditField.FontWeight = 'bold';
            app.RunTimeEditField.Tooltip = {'sec'};
            app.RunTimeEditField.Position = [915 432 100 22];

            % Create ResidenceTimeEditFieldLabel
            app.ResidenceTimeEditFieldLabel = uilabel(app.UIFigure);
            app.ResidenceTimeEditFieldLabel.HorizontalAlignment = 'right';
            app.ResidenceTimeEditFieldLabel.FontWeight = 'bold';
            app.ResidenceTimeEditFieldLabel.Position = [800 368 100 22];
            app.ResidenceTimeEditFieldLabel.Text = 'Residence  Time';

            % Create ResidenceTimeEditField
            app.ResidenceTimeEditField = uieditfield(app.UIFigure, 'numeric');
            app.ResidenceTimeEditField.FontWeight = 'bold';
            app.ResidenceTimeEditField.Tooltip = {'years'};
            app.ResidenceTimeEditField.Position = [915 368 100 22];

            % Create SequesteredCarbonLabel
            app.SequesteredCarbonLabel = uilabel(app.UIFigure);
            app.SequesteredCarbonLabel.HorizontalAlignment = 'right';
            app.SequesteredCarbonLabel.FontWeight = 'bold';
            app.SequesteredCarbonLabel.Position = [778 337 122 22];
            app.SequesteredCarbonLabel.Text = 'Sequestered Carbon';

            % Create SequesteredCarbonEditField
            app.SequesteredCarbonEditField = uieditfield(app.UIFigure, 'numeric');
            app.SequesteredCarbonEditField.FontWeight = 'bold';
            app.SequesteredCarbonEditField.Tooltip = {'PgC'};
            app.SequesteredCarbonEditField.Position = [915 337 100 22];

            % Create CarbonExportEditFieldLabel
            app.CarbonExportEditFieldLabel = uilabel(app.UIFigure);
            app.CarbonExportEditFieldLabel.HorizontalAlignment = 'right';
            app.CarbonExportEditFieldLabel.FontWeight = 'bold';
            app.CarbonExportEditFieldLabel.Position = [812 400 88 22];
            app.CarbonExportEditFieldLabel.Text = 'Carbon Export';

            % Create CarbonExportEditField
            app.CarbonExportEditField = uieditfield(app.UIFigure, 'numeric');
            app.CarbonExportEditField.FontWeight = 'bold';
            app.CarbonExportEditField.Tooltip = {'C export below the surface layer [PgC / year]'};
            app.CarbonExportEditField.Position = [915 400 100 22];

            % Create VerticalProfileDropDownLabel
            app.VerticalProfileDropDownLabel = uilabel(app.UIFigure);
            app.VerticalProfileDropDownLabel.HorizontalAlignment = 'right';
            app.VerticalProfileDropDownLabel.Position = [49 412 82 22];
            app.VerticalProfileDropDownLabel.Text = 'Vertical Profile';

            % Create VerticalProfileDropDown
            app.VerticalProfileDropDown = uidropdown(app.UIFigure);
            app.VerticalProfileDropDown.Items = {'Martin Curve', 'Vertical Migration', 'Benthos', 'Exponential'};
            app.VerticalProfileDropDown.DropDownOpeningFcn = createCallbackFcn(app, @VerticalProfileDropDownOpening, true);
            app.VerticalProfileDropDown.ValueChangedFcn = createCallbackFcn(app, @VerticalProfileDropDownValueChanged, true);
            app.VerticalProfileDropDown.ClickedFcn = createCallbackFcn(app, @VerticalProfileDropDownClicked, true);
            app.VerticalProfileDropDown.Position = [50 383 140 30];
            app.VerticalProfileDropDown.Value = 'Martin Curve';

            % Create p1EditFieldLabel
            app.p1EditFieldLabel = uilabel(app.UIFigure);
            app.p1EditFieldLabel.HorizontalAlignment = 'right';
            app.p1EditFieldLabel.Position = [223 387 25 22];
            app.p1EditFieldLabel.Text = 'p1';

            % Create p1EditField
            app.p1EditField = uieditfield(app.UIFigure, 'numeric');
            app.p1EditField.ValueChangedFcn = createCallbackFcn(app, @p1EditFieldValueChanged, true);
            app.p1EditField.Tooltip = {'Selection longitude width (in indices)'};
            app.p1EditField.Position = [254 385 47 22];
            app.p1EditField.Value = 1;

            % Create p2EditFieldLabel
            app.p2EditFieldLabel = uilabel(app.UIFigure);
            app.p2EditFieldLabel.HorizontalAlignment = 'right';
            app.p2EditFieldLabel.Position = [311 386 25 22];
            app.p2EditFieldLabel.Text = 'p2';

            % Create p2EditField
            app.p2EditField = uieditfield(app.UIFigure, 'numeric');
            app.p2EditField.ValueChangedFcn = createCallbackFcn(app, @p2EditFieldValueChanged, true);
            app.p2EditField.Tooltip = {'Selection latitude width (in indices)'};
            app.p2EditField.Position = [346 384 47 22];
            app.p2EditField.Value = 1;

            % Create AreaEditFieldLabel
            app.AreaEditFieldLabel = uilabel(app.UIFigure);
            app.AreaEditFieldLabel.HorizontalAlignment = 'right';
            app.AreaEditFieldLabel.FontWeight = 'bold';
            app.AreaEditFieldLabel.Position = [633 432 32 22];
            app.AreaEditFieldLabel.Text = 'Area';

            % Create AreaEditField
            app.AreaEditField = uieditfield(app.UIFigure, 'numeric');
            app.AreaEditField.ValueDisplayFormat = '%11g';
            app.AreaEditField.FontWeight = 'bold';
            app.AreaEditField.Position = [680 432 87 22];

            % Create secLabel
            app.secLabel = uilabel(app.UIFigure);
            app.secLabel.Position = [1023 432 25 22];
            app.secLabel.Text = 'sec';

            % Create PgCyearLabel
            app.PgCyearLabel = uilabel(app.UIFigure);
            app.PgCyearLabel.Position = [1024 400 55 22];
            app.PgCyearLabel.Text = 'PgC/year';

            % Create yearLabel
            app.yearLabel = uilabel(app.UIFigure);
            app.yearLabel.Position = [1025 368 28 22];
            app.yearLabel.Text = 'year';

            % Create PgCLabel
            app.PgCLabel = uilabel(app.UIFigure);
            app.PgCLabel.Position = [1025 337 28 22];
            app.PgCLabel.Text = 'PgC';

            % Create km2Label
            app.km2Label = uilabel(app.UIFigure);
            app.km2Label.Position = [772 432 33 22];
            app.km2Label.Text = 'km^2';

            % Create Image
            app.Image = uiimage(app.UIFigure);
            app.Image.Position = [1037 648 106 94];
            app.Image.ImageSource = fullfile(pathToMLAPP, 'Picture2.png');

            % Create gCm2yearLabel
            app.gCm2yearLabel = uilabel(app.UIFigure);
            app.gCm2yearLabel.Position = [130 441 73 22];
            app.gCm2yearLabel.Text = 'gC/m^2/year';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = TransportMatrixAppSeaQuester_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end