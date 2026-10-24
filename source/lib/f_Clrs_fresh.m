function Colors = f_Clrs_fresh(n, NW)
%Alexey Ryabov
%give n fresh distinct colors
% https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
%'#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'

arguments
    n double
    NW.Order  {mustBeMember(NW.Order,{'my','convenient', 'rainbow'})} = 'my'
end

if n < 21
    switch NW.Order
        case 'my'
            Colors = [...
                204, 51, 0;
                0, 153, 255
                245, 130, 48;...% OrangeR: 245G: 130B: 48
                145, 30,  180;...% PurpleR: 145G: 30B: 180
                60,  180, 75;...% GreenR: 60, 180, 75
                170, 110, 40;...% BrownR: 170G: 110B: 40
                255, 200, 25; %Darker Yellow %255, 225, 25;...% YellowR: 255G: 225B: 25
                230, 190, 255;...% LavenderR: 230, 190, 255
                210, 245, 60;...% LimeR: 210, 245, 60
                230, 25, 75;
                70, 240, 240;
                240, 50, 230;
                250, 190, 212;
                255, 250, 200;
                0,   128, 128;...% TealR: 0G: 128B: 128
                170, 255, 195;
                0, 0, 128;
                128, 128, 0;
                255, 215, 180;
                0,   130, 200;...% BlueR: 0G: 130B: 200
                128, 128, 128;
                0, 0, 0
                ]/255;
            Colors = Colors(1:n, :);
            
               Colors = rgb2hsv(Colors);
            [~, ind] = sort(Colors(:, 1));
            Colors = Colors(ind, :);
            Colors = hsv2rgb(Colors );
        case 'convenient'
            Colors = [...
                230, 25, 75;
                60, 180, 75;...
                255, 200, 25; %Darker Yellow %255, 225, 25;
                0, 0, 128; ...%Navi
                245, 130, 48;...% OrangeR: 245G: 130B: 48
                145, 30,  180;...% PurpleR: 145G: 30B: 180
                70, 240, 240;
                240, 50, 230; %Magenta
                210, 245, 60;...% LimeR: 210, 245, 60
                250,190,212; ... %Pink
                128, 0, 0;... %Maroon
                230, 190, 255;...% LavenderR: 230, 190, 255
                170, 110, 40;...% BrownR: 170G: 110B: 40
                255, 250, 200;... %Beige
                0,   128, 128;...% TealR: 0G: 128B: 128
                170, 255, 195;...%Mint
                128, 128, 0;...%Olive
                255, 215, 180;... %Apricot
                0,   130, 200;...% BlueR: 0G: 130B: 200
                128, 128, 128; ...%grey
                0, 0, 0
                ]/255;
            Colors = Colors(1:n, :);  
            case 'rainbow'
               Colors = [...
                230, 25, 75;
                60, 180, 75;...
                255, 200, 25; %Darker Yellow %255, 225, 25;
                0,   130, 200;...% BlueR: 0G: 130B: 200
                245, 130, 48;...% OrangeR: 245G: 130B: 48
                145, 30,  180;...% PurpleR: 145G: 30B: 180
                70, 240, 240;
                240, 50, 230; %Magenta
                210, 245, 60;...% LimeR: 210, 245, 60
                250,190,212; ... %Pink
                0,   128, 128;...% TealR: 0G: 128B: 128
                230, 190, 255;...% LavenderR: 230, 190, 255
                170, 110, 40;...% BrownR: 170G: 110B: 40
                255, 250, 200;... %Beige
                128, 0, 0;... %Maroon
                170, 255, 195;...%Mint
                128, 128, 0;...%Olive
                255, 215, 180;... %Apricot
                0, 0, 128; ...%Navi
                ]/255;
            Colors = rgb2hsv(Colors);
            [~, ind] = sort(Colors(:, 3));
            Colors = Colors(ind(1:n), :);
             
            [~, ind] = sort(Colors(:, 1));
            Colors = Colors(ind, :);
            Colors = hsv2rgb(Colors );
    end
else
    Colors = jet(n);
end
