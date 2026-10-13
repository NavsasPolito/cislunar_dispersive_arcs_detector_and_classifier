function [] = saveFigures(picFolder,pdfBool,pngBool,figBool)
% saveFigures saves all the opened figures in fig,pdf and png format
% inside the chosen directory.
%
%   EXAMPLE USAGE:
%       saveFigures(picFolder,pdfBool,pngBool,figBool)
%
%
%   AUTHOR:
%       Francesco Fiorina, NavSAS - DET - Politecnico di Torino, 14/01/2026
%
%   VERSION:
%       1.0
%
%   COPYRIGHT:
%       NavSAS - DET - Politecnico di Torino, 2026. All rights reserved.
    %dbstop if error
if nargin < 1
    picFolder = [];
    pdfBool = 1;
    pngBool = 1;
    figBool = 1;
end

if nargin == 1
    pdfBool = 1;
    pngBool = 1;
    figBool = 1;
end

if ~exist('picFolder','var') || isempty(picFolder)
    picFolder = uigetdir('../');
end

% Directory where to save
saveDir = picFolder;

% Find all open figures
figs = findobj('Type','figure');
if isempty(figs)
    return
end
figNums = [figs.Number];

for k = min(figNums):max(figNums)

    % Activate the figure (optional but clean)
    fig = figure(k);

    % Extract figure name
    figName = replace(fig.Name,{'\','$','{','}',' ','/'},'_');

    % Fallback if empty (MATLAB default)
    if isempty(figName)
        figName = sprintf('Figure%d', figs(k).Number);
    end

    % Build paths
    fileNameFig = fullfile(saveDir, figName + ".fig");
    fileNamePDF = fullfile(saveDir, figName + ".pdf");
    fileNamePNG = fullfile(saveDir, figName + ".png");

    formatFig=string([]);
    formatPDF=string([]);
    formatPng=string([]);

    %OPTIONAL: Save .fig
    if figBool 
        savefig(fig, fileNameFig);
        formatFig = 'fig';
    end

    % Save PDF (vector)
    if pdfBool
        exportgraphics(fig, fileNamePDF, ...
            'ContentType','vector', ...
            'BackgroundColor','none');
            formatPDF=' PDF ';
    end

    % Save PNG (raster)
    if pngBool
        exportgraphics(fig, fileNamePNG, ...
            'Resolution', 300, ...             % 300 dpi
            'BackgroundColor','white');        % best for PNG
            formatPng='png ';
    end

    fprintf('Saved: %s (%s%s%s)\n', figName,formatFig,formatPDF,formatPng);
end

end