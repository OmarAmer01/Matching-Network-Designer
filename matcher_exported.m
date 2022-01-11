classdef matcher_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        UITable                        matlab.ui.control.Table
        StubPropertiesPanel            matlab.ui.container.Panel
        sR_0                           matlab.ui.control.NumericEditField
        CharachteristicImpedanceLabel_3  matlab.ui.control.Label
        StubTerminationImpedanceButtonGroup  matlab.ui.container.ButtonGroup
        slabel                         matlab.ui.control.Label
        ShortCircuitTerminationButton  matlab.ui.control.RadioButton
        OpenCircuitTerminationButton   matlab.ui.control.RadioButton
        TransmissionLinePropertiesPanel_2  matlab.ui.container.Panel
        tlR0                           matlab.ui.control.NumericEditField
        CharachteristicImpedanceLabel_2  matlab.ui.control.Label
        DesignNetworkButton            matlab.ui.control.Button
        SourcePropertiesPanel          matlab.ui.container.Panel
        wavelength                     matlab.ui.control.NumericEditField
        WavelengthmLabel               matlab.ui.control.Label
        MatchingNetworkTypeDropDownLabel  matlab.ui.control.Label
        MatchingNetworkTypeDropDown    matlab.ui.control.DropDown
        LoadPropertiesPanel            matlab.ui.container.Panel
        ZrealjImagLabel                matlab.ui.control.Label
        zlImag                         matlab.ui.control.NumericEditField
        ImpedanceImaginaryPartEditFieldLabel  matlab.ui.control.Label
        zlReal                         matlab.ui.control.NumericEditField
        ImpedanceRealPartEditFieldLabel  matlab.ui.control.Label
        UIAxes                         matlab.ui.control.UIAxes
    end


    methods (Access = private)

        function [gamma, circle] = gammaCircle(~, Z_L, Z_0, impOrAdm)
            % Returns the reflection coeffecient and circle of constant
            % reflection coeffecient magnitude.

            % Normalize the load impedance with respect to the
            % charachteristic impedance of the line

            mag = abs(Z_0);
            if impOrAdm == 1
                % Impedance

                Zn_0 = Z_0 / mag;
                Zn_L = Z_L / mag;
            else
                % Admittance
                Zn_0 = (Z_0 / mag);
                Zn_L = (Z_L / mag);

            end


            % Calculate the reflection coeffecient
            gamma = ((Zn_L - Zn_0)/(Zn_L + Zn_0));
            gammaMag = abs(gamma);

            % Get the circle of constant reflection coeffecient magnitude.
            syms x y
            circle = x^2 + y^2 == gammaMag^2;
        end

        function [x, y] = enterChart(app, Z_L, R_0, impOrAdm)
            % Gets the point that represents the load impedance on the
            % Smith chart.
            if Z_L ~= inf+1i*inf
                [gamma, ~] = gammaCircle(app, Z_L, R_0, impOrAdm);
                gammaMag = abs(gamma);
                x = gammaMag * cos(angle(gamma));
                y = gammaMag * sin(angle(gamma));
            else
                x = 1;
                y = 0;
            end


        end

        function lineEqn = getStraightLine(~, a, b)
            % Gets the equation of the straight line that goes from a => b.
            m = (a(2)-b(2)) / (a(1)-b(1));
            c = b(2) - m*b(1);
            syms x y
            lineEqn = y == m*x + c;
        end

        function [Y_lx, Y_ly] = getAdmittancePt(~, Z_lx, Z_ly)
            % Gets the admittance point on the Smith chart from the point
            % of load impedance.

            Y_lx = -Z_lx;
            Y_ly = -Z_ly;

        end

        function [toX, toY, d] = moveToR1onChart(app, Z_L, R_0, impOrAdm)
            % Move on the circle of constant gamma magnitude towards the
            % generator until we intersect with the circle of r = 1.
            % Returns the points of intersection as well as the normalized
            % distance between the load point and the intersection points.

            syms x y
            constRCirc = (x-0.5)^2 + y^2 == 0.5^2;
            [~, gammaCirc] = gammaCircle(app, Z_L, R_0, impOrAdm);
            [toX, toY] = solve(constRCirc, gammaCirc);

            ResAngles = angle(toX + 1j*toY); % angle(gamma) of each point of intersection
            [xz, yz] = enterChart(app, Z_L, R_0, impOrAdm);

            if impOrAdm == 1
                % Impedance
                zlAngle = angle(xz + 1j*yz);
            else
                % Admittance
                zlAngle = angle(-xz - 1j*yz);
            end


            radialDist = zeros(2);
            for i = 1:2
                if zlAngle < ResAngles(i)
                    radialDist(i) = (2*pi)-(double(ResAngles(i)) - zlAngle);
                else
                    radialDist(i) = double(zlAngle - double(ResAngles(i)));
                end
            end



            d = (0.5*rad2deg(radialDist)/360);

        end

        function [toPlot,reac] = getNReacFromPoint(~,x, y)
            % Gets the normalized reactance of a given point on the
            % Smith chart.

            gamma = eval(x+1j*y);
            gr = real(gamma);
            gi = imag(gamma);

            syms re
            % Equation of the circle of constant reactance.
            reacCircleEqn = (gr-1).^2 + (gi-(1/re)).^2 == (1/re).^2;

            % We have two solutions
            reac(1) = solve(reacCircleEqn(1));
            reac(2) = solve(reacCircleEqn(2));

            reac = eval(reac);
            syms x y
            toPlot(1) = (x-1)^2 + (y-(1/reac(1)))^2 == (1/reac(1))^2;
            toPlot(2) = (x-1)^2 + (y-(1/reac(2)))^2 == (1/reac(2))^2;
        end

        function [toPlot,newReac] = changeLineAndGetReac(app, TLR_0, sR_0, x, y, impOrAdm)
            % Gets the normalized reactance of a given point on the
            % Smith chart after changing the line.

            [~, oldReac] = getNReacFromPoint(app, x, y);

            if impOrAdm == 1
                % Impedance
                newReac = oldReac * TLR_0 / sR_0;
            else
                newReac = oldReac * (TLR_0 / sR_0)^-1;
            end



            syms x y
            toPlot(1) = (x-1)^2 + (y-(1/newReac(1)))^2 == (1/newReac(1))^2;
            toPlot(2) = (x-1)^2 + (y-(1/newReac(2)))^2 == (1/newReac(2))^2;

        end



        function dist = getNDist(~, fromX, fromY, toX, toY, genOrLoad)
            fromAngle = angle(fromX+1j*fromY);
            toAngle = angle(toX+1j*toY);

            radialDist = zeros(2);
            for i=1:2
                if toAngle < fromAngle(i)
                    radialDist(i) = (2*pi)-(fromAngle(i) - toAngle);
                else
                    radialDist(i) = toAngle - fromAngle(i);
                end

            end


            dist = (0.5*rad2deg(radialDist)/360);
            if(genOrLoad == 0)
                % WTL scale

                dist = 0.5 - dist;

            end

        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.UITable.Visible = 'off';
            % Update the label value to match the user input
            if app.zlImag.Value < 0
                app.ZrealjImagLabel.Text = "Z = " +num2str(app.zlReal.Value) + "-j" + num2str(abs(app.zlImag.Value)) + " Ω";
            elseif app.zlImag.Value == 0
                app.ZrealjImagLabel.Text = "Z = " +num2str(app.zlReal.Value) + " Ω";

            else
                app.ZrealjImagLabel.Text = "Z = " +num2str(app.zlReal.Value) + "+j" + num2str(app.zlImag.Value) + " Ω" ;
            end
            if app.zlReal.Value == 0 && app.zlImag.Value ~= 0
                if app.zlImag.Value == 0
                    app.ZrealjImagLabel.Text = "Z = 0 Ω" ;
                elseif app.zlImag.Value < 0
                    app.ZrealjImagLabel.Text = "Z = -j" + num2str(abs(app.zlImag.Value)) + " Ω" ;
                else
                    app.ZrealjImagLabel.Text = "Z = j" + num2str(app.zlImag.Value) + " Ω" ;
                end
            end

            % Update the label value to match the user input



            app.slabel.Text = "Z = 0 Ω";

        end

        % Value changed function: zlReal
        function zlRealValueChanged(app, event)
            % Update the label value to match the user input
            if app.zlImag.Value < 0
                app.ZrealjImagLabel.Text = "Z = " +num2str(app.zlReal.Value) + "-j" + num2str(abs(app.zlImag.Value)) + " Ω";
            elseif app.zlImag.Value == 0
                app.ZrealjImagLabel.Text = "Z = " +num2str(app.zlReal.Value) + " Ω";

            else
                app.ZrealjImagLabel.Text = "Z = " +num2str(app.zlReal.Value) + "+j" + num2str(app.zlImag.Value) + " Ω" ;
            end
            if app.zlReal.Value == 0 && app.zlImag.Value ~= 0
                if app.zlImag.Value == 0
                    app.ZrealjImagLabel.Text = "Z = 0 Ω" ;
                elseif app.zlImag.Value < 0
                    app.ZrealjImagLabel.Text = "Z = -j" + num2str(abs(app.zlImag.Value)) + " Ω" ;
                else
                    app.ZrealjImagLabel.Text = "Z = j" + num2str(app.zlImag.Value) + " Ω" ;
                end
            end
        end

        % Value changed function: zlImag
        function zlImagValueChanged(app, event)
            % Update the label value to match the user input
            if app.zlImag.Value < 0
                app.ZrealjImagLabel.Text = "Z = " +num2str(app.zlReal.Value) + "-j" + num2str(abs(app.zlImag.Value)) + " Ω";
            elseif app.zlImag.Value == 0
                app.ZrealjImagLabel.Text = "Z = " +num2str(app.zlReal.Value) + " Ω";

            else
                app.ZrealjImagLabel.Text = "Z = " +num2str(app.zlReal.Value) + "+j" + num2str(app.zlImag.Value) + " Ω" ;
            end
            if app.zlReal.Value == 0 && app.zlImag.Value ~= 0
                if app.zlImag.Value == 0
                    app.ZrealjImagLabel.Text = "Z = 0 Ω" ;
                elseif app.zlImag.Value < 0
                    app.ZrealjImagLabel.Text = "Z = -j" + num2str(abs(app.zlImag.Value)) + " Ω" ;
                else
                    app.ZrealjImagLabel.Text = "Z = j" + num2str(app.zlImag.Value) + " Ω" ;
                end
            end
        end

        % Button pushed function: DesignNetworkButton
        function DesignNetworkButtonPushed(app, event)
            % Clear Plot
            cla(app.UIAxes)

            % Design the matching network.

            % Acquire the parameters

            Z_L = app.zlReal.Value + 1j * app.zlImag.Value; % Load Impedance
            TLR_0 = app.tlR0.Value; % Transmission Line Charachteristic Impedance
            stubR_0 = app.sR_0.Value;

            lambda = app.wavelength.Value;
            userSel = app.MatchingNetworkTypeDropDown.Value;


            switch userSel
                case "Single Stub - Series"
                    % Aquire the circle of constant gamma magnitude and
                    % plot it.
                    [~, circleEqn] = gammaCircle(app, Z_L, TLR_0, 1);
                    fimplicit(app.UIAxes, circleEqn, 'Color',[0 0 0])
                    hold(app.UIAxes, 'on')

                    % Enter The Smith chart from the Z_L point and plot it.
                    [Z_lx, Z_ly] = enterChart(app, Z_L, TLR_0,1);
                    plot(app.UIAxes,Z_lx, Z_ly,'-o', "Color",[0 0 1], "LineWidth",3)


                    % Find and plot the points of intersection of the
                    % circle of constant gamma magnitude and the
                    % circle of r = 1.
                    % Also returns the normalized disance moved towards
                    % the generator.
                    [interX, interY, d] = moveToR1onChart(app, Z_L, TLR_0,1);

                    % To de-normalize
                    d = d*lambda;
                    plot(app.UIAxes,interX, interY, 'X', 'LineWidth',3)

                    % Find the reactance of the points of intersection.
                    [reacCircle, ~] = changeLineAndGetReac(app, TLR_0, stubR_0, interX, interY, 1);
                    fimplicit(app.UIAxes, reacCircle)



                    if app.ShortCircuitTerminationButton.Value == true
                        % Enter the chart from the short circuit point
                        [~, scCircle] = gammaCircle(app, 0, TLR_0,1);
                        [x_c, y_c] = enterChart(app, 0, TLR_0,1);

                        % Plot the SC circle and the SC point
                        fimplicit(app.UIAxes,scCircle, 'Color',[0 0 0])
                        plot(app.UIAxes, x_c, y_c, 'X', 'LineWidth',3,"Color",[1 0 1])
                    else
                        % Enter the chart from the open circuit point
                        [~, scCircle] = gammaCircle(app, 0, TLR_0,1);
                        [x_c, y_c] = enterChart(app, (inf+inf*1i), TLR_0,1);

                        % Plot the SC circle and the SC point
                        fimplicit(app.UIAxes,scCircle, 'Color',[0 0 0])
                        plot(app.UIAxes, x_c, y_c, 'X', 'LineWidth',3,"Color",[1 0 1])
                    end


                    % Find the points of intersection of the const
                    % reactance circle and the SC circle as well as the
                    % normalized distance moved towards the generator.
                    % We also plot the points.

                    % Refuse the first solution
                    [ansX1, ansY1] = solve(reacCircle(2), scCircle);
                    stubLenX(1) = eval(ansX1(2));
                    stubLenY(1) = eval(ansY1(2));

                    [ansX2, ansY2] = solve(reacCircle(1), scCircle);
                    stubLenX(2) = eval(ansX2(2));
                    stubLenY(2) = eval(ansY2(2));

                    % Plot the points
                    plot(app.UIAxes, stubLenX, stubLenY, 'X', 'LineWidth',3,"Color",[0.5 1 0.5])
                    %plot(app.UIAxes, stubLenX2, stubLenY2, 'X', 'LineWidth',3,"Color",[0.5 1 0.5])


                    stubLen = lambda * getNDist(app,stubLenX,stubLenY,x_c,y_c,1);

                    % Display The Results.
                    app.UITable.Data = [[d(1) stubLen(1)]; [d(2) stubLen(2)]];
                    app.UITable.Visible = 'on';


                case "Single Stub - Parallel"
                    % Aquire the circle of constant gamma magnitude and
                    % plot it.
                    [~, circleEqn] = gammaCircle(app, Z_L, TLR_0,0);
                    fimplicit(app.UIAxes, circleEqn, 'Color',[0 0 0])
                    hold(app.UIAxes, 'on')

                    % Enter The Smith chart from the Z_L point and plot it.
                    [Z_lx, Z_ly] = enterChart(app, Z_L, TLR_0,0);
                    plot(app.UIAxes,Z_lx, Z_ly,'-o', "Color",[0 0 1], "LineWidth",3)


                    %Get the admittance point and plot it.
                    [Y_lx, Y_ly] = getAdmittancePt(app, Z_lx, Z_ly);
                    plot(app.UIAxes, Y_lx, Y_ly, '-o', "Color",[1 0 0], "LineWidth",3)

                    % Find and plot the points of intersection of the
                    % circle of constant gamma magnitude and the
                    % circle of g = 1.
                    % Also returns the normalized disance moved towards
                    % the generator.

                    %Y_L = 1/Z_L;
                    [interX, interY, d] = moveToR1onChart(app, Z_L, TLR_0,0);

                    % To de-normalize
                    d = d*lambda;
                    plot(app.UIAxes,interX, interY, 'X', 'LineWidth',3)

                    % Find the reactance of the points of intersection.
                    [reacCircle, ~] = changeLineAndGetReac(app, TLR_0, stubR_0, interX, interY, 0);
                    fimplicit(app.UIAxes, reacCircle)



                    if app.ShortCircuitTerminationButton.Value == true
                        % Enter the chart from the open circuit point
                        [~, scCircle] = gammaCircle(app, 0, TLR_0,0);
                        [x_c, y_c] = enterChart(app, (inf+1i*inf), TLR_0,0);

                        % Plot the SC circle and the SC point
                        fimplicit(app.UIAxes,scCircle, 'Color',[0 0 0])
                        plot(app.UIAxes, x_c, y_c, 'X', 'LineWidth',3,"Color",[1 0 1])
                    else
                        % Enter the chart from the short circuit point
                        [~, scCircle] = gammaCircle(app, 0, TLR_0,0);
                        [x_c, y_c] = enterChart(app, 0, TLR_0,0);

                        % Plot the SC circle and the SC point
                        fimplicit(app.UIAxes,scCircle, 'Color',[0 0 0])
                        plot(app.UIAxes, x_c, y_c, 'X', 'LineWidth',3,"Color",[1 0 1])
                    end


                    % Find the points of intersection of the const
                    % reactance circle and the SC circle as well as the
                    % normalized distance moved towards the generator.
                    % We also plot the points.

                    % Refuse the first solution
                    [ansX1, ansY1] = solve(reacCircle(2), scCircle);
                    stubLenX(1) = eval(ansX1(2));
                    stubLenY(1) = eval(ansY1(2));

                    [ansX2, ansY2] = solve(reacCircle(1), scCircle);
                    stubLenX(2) = eval(ansX2(2));
                    stubLenY(2) = eval(ansY2(2));

                    % Plot the points
                    plot(app.UIAxes, stubLenX, stubLenY, 'X', 'LineWidth',3,"Color",[0.5 1 0.5])
                    %plot(app.UIAxes, stubLenX2, stubLenY2, 'X', 'LineWidth',3,"Color",[0.5 1 0.5])


                    stubLen = lambda * getNDist(app,stubLenX,stubLenY,x_c,y_c,1);

                    % Display The Results.
                    app.UITable.Data = [[d(1) stubLen(1)]; [d(2) stubLen(2)]];
                    app.UITable.Visible = 'on';

                case "Double Stub - Parallel"
                    stubDist = app.sR_0.Value;




                    % Enter The Smith chart from the Z_L point and plot it.
                    [Z_lx, Z_ly] = enterChart(app, Z_L, TLR_0,0);
                    %plot(app.UIAxes,Z_lx, Z_ly,'-o', "Color",[0 0 0], "LineWidth",3)
                    hold(app.UIAxes, 'on')


                    [~, scCircle] = gammaCircle(app, 0, TLR_0,1);

                    % Plot the SC circle
                    fimplicit(app.UIAxes,scCircle, 'Color',[0 0 0])

                    %Get the admittance point and plot it.
                    [Y_lx, Y_ly] = getAdmittancePt(app, Z_lx, Z_ly);
                    plot(app.UIAxes, Y_lx, Y_ly, '-o', "Color",[1 0 0], "LineWidth",3)

                    syms x y u v a


                    angleOfRot = deg2rad(-stubDist*90/0.125);
                    rotMat = [[cos(a) -sin(a)];[sin(a) cos(a)]];
                    newAxes = rotMat * [x ; y];
                    u = newAxes(1);
                    v = newAxes(2);

                    G1 = (x-0.5)^2 + y^2 == 0.5^2;
                    fimplicit(app.UIAxes,subs(G1, a, angleOfRot))
                    rotatedG1 = (u-0.5)^2 + v^2 == 0.5^2;
                    fimplicit(app.UIAxes,subs(rotatedG1, a, angleOfRot))
                    zn = (app.zlReal.Value+j*app.zlImag.Value) / app.tlR0.Value;
                    r = real(zn^-1);

                    circleOfConstG = (x - (r/(1+r)))^2 + y^2 == (1/(1+r))^2;
                    fimplicit(app.UIAxes,circleOfConstG)
                    [x_ya, y_ya] = solve(circleOfConstG, rotatedG1);
                    x_ya = (subs(x_ya, a, angleOfRot));
                    y_ya = (subs(y_ya, a, angleOfRot));



                    plot(app.UIAxes, x_ya, y_ya,'X' ,"Color",[0 1 0] ,"LineWidth",3)

                    [reacArc, reac] = getNReacFromPoint(app, x_ya, y_ya);
                    %fimplicit(app.UIAxes, reacArc,"LineStyle","--")

                    stubSus = reac - imag(zn^-1);
                    susArc(1) = (x-1)^2 + (y-(1/stubSus(1)))^2 == (1/stubSus(1))^2;
                    susArc(2) = (x-1)^2 + (y-(1/stubSus(2)))^2 == (1/stubSus(2))^2;

                    %fimplicit(app.UIAxes, susArc(1))
                    %fimplicit(app.UIAxes, susArc(2))

                    % now get the suseptance on big gamma circle
                    [susX1tmp, susY1tmp] = solve(scCircle, susArc(1));
                    [susX2tmp, susY2tmp] = solve(scCircle, susArc(2));

                    susX1 = eval(susX1tmp);
                    susX2 = eval(susX2tmp);
                    susY1 = eval(susY1tmp);
                    susY2 = eval(susY2tmp);

                    % plot(app.UIAxes, susX1(2), susY1(2),'X' ,"Color",[0 0 1] ,"LineWidth",3)
                    % plot(app.UIAxes, susX2(2), susY2(2),'X' ,"Color",[0 0 1] ,"LineWidth",3)

                    xsc = 1;
                    ysc = 0;

                    dist1 = getNDist(app,susX2, susY2, xsc, ysc, 1);
                    dist2 = getNDist(app,susX1, susY1, xsc, ysc, 1);
                    app.UITable.ColumnName = {'l_1'; 'l_2'};
                    app.UITable.Visible = 'On';




                    % we use this to calc the 2nd solution

                    gammaCirc = x^2 + y^2 == x_ya(1).^2 + y_ya(1).^2;
                    
                    fimplicit(app.UIAxes,gammaCirc,'DisplayName','GAMMA CIRCLE','Color',[1 0 1]);

                    [x_yb, y_yb] = solve(gammaCirc, G1);
                    plot(app.UIAxes, x_yb, y_yb, 'o')

                    [~, susBtmp] = getNReacFromPoint(app, x_yb, y_yb);
                    susB = susBtmp- imag(zn^-1) ;
                    susBArc(1) = (x-1)^2 + (y-(1/susB(1)))^2 == (1/susB(1))^2;
                    susBArc(2) = (x-1)^2 + (y-(1/susB(2)))^2 == (1/susB(2))^2;
                    
                    gammaCircB = x^2 + y^2 == x_yb(1).^2 + y_yb(1).^2;

                    
                    
                    [x_yb2, y_yb2] = solve(gammaCircB, G1);
                    plot(app.UIAxes, x_yb2, y_yb2, 'o')

                    [~, susBtmp2] = getNReacFromPoint(app, x_yb, y_yb);
                    susB2 = susBtmp2- imag(zn^-1) ;
                    susBArc2(1) = (x-1)^2 + (y-(1/susB2(1)))^2 == (1/susB2(1))^2;
                    susBArc2(2) = (x-1)^2 + (y-(1/susB2(2)))^2 == (1/susB2(2))^2;
                    
                    
                    % now get the suseptance on big gamma circle
                    [susX1tmp2, susY1tmp2] = solve(scCircle, susBArc(1)); %%%%%%
                    [susX2tmp2, susY2tmp2] = solve(scCircle, susBArc(2));

                    sus2X1 = eval(susX1tmp2);
                    sus2X2 = eval(susX2tmp2);
                    sus2Y1 = eval(susY1tmp2);
                    sus2Y2 = eval(susY2tmp2);

                    plot(app.UIAxes, sus2X1(1), sus2Y1(1),'*' ,"Color",[1 0 0] ,"LineWidth",3)
                    plot(app.UIAxes, sus2X2(2), sus2Y2(2),'*' ,"Color",[1 0 0] ,"LineWidth",3)

                    dist11 = getNDist(app,sus2X2, sus2Y2, xsc, ysc, 1);
                    dist22 = getNDist(app,sus2X1, sus2Y1, xsc, ysc, 1);


                    app.UITable.Data = [[dist2(2) dist11(2)]; [dist1(2) dist22(1)]];
            end



        end

        % Callback function
        function sRValueChanged(app, event)
            % Update the label value to match the user input
            if app.sI.Value < 0
                app.slabel.Text = "Z = " +num2str(app.sR.Value) + "-j" + num2str(abs(app.sI.Value)) + " Ω";
            elseif app.sI.Value == 0
                app.slabel.Text = "Z = " +num2str(app.sR.Value) + " Ω";

            else
                app.slabel.Text = "Z = " +num2str(app.sR.Value) + "+j" + num2str(app.sI.Value) + " Ω" ;
            end
            if app.sR.Value == 0 && app.sI.Value ~= 0
                if app.sI.Value == 0
                    app.slabel.Text = "Z = 0 Ω" ;
                elseif app.sI.Value < 0
                    app.slabel.Text = "Z = -j" + num2str(abs(app.sI.Value)) + " Ω" ;
                else
                    app.slabel.Text = "Z = j" + num2str(app.sI.Value) + " Ω" ;
                end
            end
        end

        % Callback function
        function sIValueChanged(app, event)
            % Update the label value to match the user input
            if app.sI.Value < 0
                app.slabel.Text = "Z = " +num2str(app.sR.Value) + "-j" + num2str(abs(app.sI.Value)) + " Ω";
            elseif app.sI.Value == 0
                app.slabel.Text = "Z = " +num2str(app.sR.Value) + " Ω";

            else
                app.slabel.Text = "Z = " +num2str(app.sR.Value) + "+j" + num2str(app.sI.Value) + " Ω" ;
            end
            if app.sR.Value == 0 && app.sI.Value ~= 0
                if app.sI.Value == 0
                    app.slabel.Text = "Z = 0 Ω" ;
                elseif app.sI.Value < 0
                    app.slabel.Text = "Z = -j" + num2str(abs(app.sI.Value)) + " Ω" ;
                else
                    app.slabel.Text = "Z = j" + num2str(app.sI.Value) + " Ω" ;
                end
            end
        end

        % Selection changed function: 
        % StubTerminationImpedanceButtonGroup
        function StubTerminationImpedanceButtonGroupSelectionChanged(app, event)
            selectedButton = app.StubTerminationImpedanceButtonGroup.SelectedObject;
            if selectedButton == app.ShortCircuitTerminationButton

                app.slabel.Text = "Z = 0 Ω";

            elseif selectedButton == app.OpenCircuitTerminationButton

                app.slabel.Text = "Z = ∞ Ω";

            end
        end

        % Value changed function: MatchingNetworkTypeDropDown
        function MatchingNetworkTypeDropDownValueChanged(app, event)
            value = app.MatchingNetworkTypeDropDown.Value;
            if value == "Double Stub - Parallel"
                app.CharachteristicImpedanceLabel_3.Text = "Normalized stub Distance: ";
                app.ShortCircuitTerminationButton.Value = true;
                app.StubTerminationImpedanceButtonGroup.Enable = 'off';
            else
                app.CharachteristicImpedanceLabel_3.Text = "Charachteristic Impedance (Ω)";
                              
                app.StubTerminationImpedanceButtonGroup.Enable = 'on';
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1138 446];
            app.UIFigure.Name = 'MATLAB App';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, '\Gamma_r')
            ylabel(app.UIAxes, '\Gamma_i')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.XLim = [-1 1];
            app.UIAxes.YLim = [-1 1];
            app.UIAxes.XAxisLocation = 'origin';
            app.UIAxes.XTick = [];
            app.UIAxes.YAxisLocation = 'origin';
            app.UIAxes.YTick = [];
            app.UIAxes.XGrid = 'on';
            app.UIAxes.XMinorGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.YMinorGrid = 'on';
            app.UIAxes.Position = [738 24 369 380];

            % Create LoadPropertiesPanel
            app.LoadPropertiesPanel = uipanel(app.UIFigure);
            app.LoadPropertiesPanel.Title = 'Load Properties';
            app.LoadPropertiesPanel.Position = [42 256 323 148];

            % Create ImpedanceRealPartEditFieldLabel
            app.ImpedanceRealPartEditFieldLabel = uilabel(app.LoadPropertiesPanel);
            app.ImpedanceRealPartEditFieldLabel.HorizontalAlignment = 'right';
            app.ImpedanceRealPartEditFieldLabel.Position = [12 96 118 22];
            app.ImpedanceRealPartEditFieldLabel.Text = 'Impedance Real Part';

            % Create zlReal
            app.zlReal = uieditfield(app.LoadPropertiesPanel, 'numeric');
            app.zlReal.ValueDisplayFormat = '%11.9g';
            app.zlReal.ValueChangedFcn = createCallbackFcn(app, @zlRealValueChanged, true);
            app.zlReal.Position = [173 96 144 22];
            app.zlReal.Value = 30;

            % Create ImpedanceImaginaryPartEditFieldLabel
            app.ImpedanceImaginaryPartEditFieldLabel = uilabel(app.LoadPropertiesPanel);
            app.ImpedanceImaginaryPartEditFieldLabel.HorizontalAlignment = 'right';
            app.ImpedanceImaginaryPartEditFieldLabel.Position = [12 55 146 22];
            app.ImpedanceImaginaryPartEditFieldLabel.Text = 'Impedance Imaginary Part';

            % Create zlImag
            app.zlImag = uieditfield(app.LoadPropertiesPanel, 'numeric');
            app.zlImag.ValueDisplayFormat = '%11.9g';
            app.zlImag.ValueChangedFcn = createCallbackFcn(app, @zlImagValueChanged, true);
            app.zlImag.Position = [173 55 145 22];
            app.zlImag.Value = 40;

            % Create ZrealjImagLabel
            app.ZrealjImagLabel = uilabel(app.LoadPropertiesPanel);
            app.ZrealjImagLabel.Position = [15 9 303 22];
            app.ZrealjImagLabel.Text = 'Z =  real+jImag';

            % Create MatchingNetworkTypeDropDown
            app.MatchingNetworkTypeDropDown = uidropdown(app.UIFigure);
            app.MatchingNetworkTypeDropDown.Items = {'Single Stub - Series', 'Single Stub - Parallel', 'Double Stub - Parallel', ''};
            app.MatchingNetworkTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @MatchingNetworkTypeDropDownValueChanged, true);
            app.MatchingNetworkTypeDropDown.Position = [538 214 178 22];
            app.MatchingNetworkTypeDropDown.Value = 'Single Stub - Series';

            % Create MatchingNetworkTypeDropDownLabel
            app.MatchingNetworkTypeDropDownLabel = uilabel(app.UIFigure);
            app.MatchingNetworkTypeDropDownLabel.HorizontalAlignment = 'right';
            app.MatchingNetworkTypeDropDownLabel.Position = [392 214 131 22];
            app.MatchingNetworkTypeDropDownLabel.Text = 'Matching Network Type';

            % Create SourcePropertiesPanel
            app.SourcePropertiesPanel = uipanel(app.UIFigure);
            app.SourcePropertiesPanel.Title = 'Source Properties';
            app.SourcePropertiesPanel.Position = [392 256 323 70];

            % Create WavelengthmLabel
            app.WavelengthmLabel = uilabel(app.SourcePropertiesPanel);
            app.WavelengthmLabel.HorizontalAlignment = 'right';
            app.WavelengthmLabel.Position = [6 9 90 22];
            app.WavelengthmLabel.Text = 'Wavelength (m)';

            % Create wavelength
            app.wavelength = uieditfield(app.SourcePropertiesPanel, 'numeric');
            app.wavelength.ValueDisplayFormat = '%11.9g';
            app.wavelength.Position = [216 9 100 22];
            app.wavelength.Value = 1;

            % Create DesignNetworkButton
            app.DesignNetworkButton = uibutton(app.UIFigure, 'push');
            app.DesignNetworkButton.ButtonPushedFcn = createCallbackFcn(app, @DesignNetworkButtonPushed, true);
            app.DesignNetworkButton.Position = [392 166 323 22];
            app.DesignNetworkButton.Text = 'Design Network';

            % Create TransmissionLinePropertiesPanel_2
            app.TransmissionLinePropertiesPanel_2 = uipanel(app.UIFigure);
            app.TransmissionLinePropertiesPanel_2.Title = 'Transmission Line Properties';
            app.TransmissionLinePropertiesPanel_2.Position = [393 334 323 70];

            % Create CharachteristicImpedanceLabel_2
            app.CharachteristicImpedanceLabel_2 = uilabel(app.TransmissionLinePropertiesPanel_2);
            app.CharachteristicImpedanceLabel_2.HorizontalAlignment = 'right';
            app.CharachteristicImpedanceLabel_2.Position = [5 9 168 22];
            app.CharachteristicImpedanceLabel_2.Text = 'Charachteristic Impedance (Ω)';

            % Create tlR0
            app.tlR0 = uieditfield(app.TransmissionLinePropertiesPanel_2, 'numeric');
            app.tlR0.ValueDisplayFormat = '%11.9g';
            app.tlR0.Position = [216 9 100 22];
            app.tlR0.Value = 75;

            % Create StubPropertiesPanel
            app.StubPropertiesPanel = uipanel(app.UIFigure);
            app.StubPropertiesPanel.Title = 'Stub Properties';
            app.StubPropertiesPanel.Position = [42 61 323 175];

            % Create StubTerminationImpedanceButtonGroup
            app.StubTerminationImpedanceButtonGroup = uibuttongroup(app.StubPropertiesPanel);
            app.StubTerminationImpedanceButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @StubTerminationImpedanceButtonGroupSelectionChanged, true);
            app.StubTerminationImpedanceButtonGroup.Title = 'Stub Termination Impedance';
            app.StubTerminationImpedanceButtonGroup.Position = [9 10 305 96];

            % Create OpenCircuitTerminationButton
            app.OpenCircuitTerminationButton = uiradiobutton(app.StubTerminationImpedanceButtonGroup);
            app.OpenCircuitTerminationButton.Text = 'Open Circuit Termination';
            app.OpenCircuitTerminationButton.Position = [8 47 154 22];

            % Create ShortCircuitTerminationButton
            app.ShortCircuitTerminationButton = uiradiobutton(app.StubTerminationImpedanceButtonGroup);
            app.ShortCircuitTerminationButton.Text = 'Short Circuit Termination';
            app.ShortCircuitTerminationButton.Position = [8 26 153 22];
            app.ShortCircuitTerminationButton.Value = true;

            % Create slabel
            app.slabel = uilabel(app.StubTerminationImpedanceButtonGroup);
            app.slabel.Position = [8 5 301 22];
            app.slabel.Text = 'Z =  real+jImag';

            % Create CharachteristicImpedanceLabel_3
            app.CharachteristicImpedanceLabel_3 = uilabel(app.StubPropertiesPanel);
            app.CharachteristicImpedanceLabel_3.HorizontalAlignment = 'right';
            app.CharachteristicImpedanceLabel_3.Position = [6 126 168 22];
            app.CharachteristicImpedanceLabel_3.Text = 'Charachteristic Impedance (Ω)';

            % Create sR_0
            app.sR_0 = uieditfield(app.StubPropertiesPanel, 'numeric');
            app.sR_0.ValueDisplayFormat = '%11.9g';
            app.sR_0.Position = [217 126 100 22];
            app.sR_0.Value = 50;

            % Create UITable
            app.UITable = uitable(app.UIFigure);
            app.UITable.BackgroundColor = [1 1 1];
            app.UITable.ColumnName = {'Stub distance from load (m)'; 'Stub length (m)'};
            app.UITable.RowName = {'Point 1,'; 'Point 2'};
            app.UITable.Position = [392 61 323 89];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = matcher_exported

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