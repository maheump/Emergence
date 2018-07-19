function [ xbeg, ybeg, xend, yend ] = DisplayCenteredText( windowPtr, text, crossX, crossY, col )

[w, h] = RectSize(Screen('TextBounds', windowPtr, text));
xbeg = round(crossX-w/2);
ybeg = round(crossY-h/2);
xend = round(crossX+w/2);
yend = round(crossY+h/2);
Screen('DrawText', windowPtr, text, xbeg, ybeg, col);

end