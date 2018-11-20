function figtopdf(myfigname, nocomment)
%FIGTOPDF - Save the current figure as a .pdf image with the (optional)
%given string (without the '.pdf' extension). If no name is given as
%argument, the current figure is saved with the default name 'myfigure'. A
%boolean can also be passed as an argument so the function doesn't display any
%comment in the command window. If multiple figures are displayed in the 
%Matlab environnement, only the last active figure that will be saved as pdf.
%
% If multiple figures are displayed in the Matlab environnement, you can
% change the figure to be saved by clicking on the figure, so the
% figure become the current figure.
%
% Syntax: figtopdf(myfigname, nocomment)
%
% Inputs:
%    myfigname - Optionnal. Name of the image to save (witouth the '.pdf' extension)
%    nocomment - optinnal. Boolean that is true if you don't want the
%                function to display any text on the cmmand window. Default is false
%               (the function will display a comment in the command window by default)
%
% Outputs:
%    none
%
% Example: 
%    figtopdf()                       % Save the current figure as "myfigure.pdf" and display a confirmation message in the comand window.
%    figtopdf('nameofmyfigure')       % Save the current figure as "nameofmyfigure.pdf" and display a confirmation message in the comand window.
%    figtopdf('nameofmyfigure', true) % Save the current figure as "nameofmyfigure.pdf" and don't display any message in the cvommand window.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Alexandre Willame & Antoine Berthelemot
% April 2013; Last revision: 29-April-2013

%------------- BEGIN CODE --------------

if (nargin == 0)                   % Multiple argument handling
    filename = 'myfigure';
else
    filename = myfigname;
end
if (nargin ~= 2)
    nocomment = false;
end

print('-dpdf','-tiff', filename);  % fig to eps
%system(['epstopdf ' filename '.eps']); % eps to pdf
%delete([filename, '.eps']);         % delete extra '.eps' file

if (~nocomment)
    str = sprintf('Image "%s.pdf" saved.', filename);
    disp(str);                      % Confirmation message in the command window
end

%------------- END OF CODE --------------
end

