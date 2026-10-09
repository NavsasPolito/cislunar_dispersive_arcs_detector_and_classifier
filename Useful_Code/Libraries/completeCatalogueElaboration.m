function  completeCatalogueElaboration(jj,settings)

if jj ==20
    % --- Command Window Log
    fprintf('-> Complete catalogue elaboration...\n')
    close all

    % Dopo aver aggiornato exportCatalogues e rilanciato l'aggregazione:
    T = exportCatalogues(settings.resultPath,[],settings.opTablePath);
    disp('Saved: complete catalogue')
    plotCataloguesGrazingAltitudeStats(T);
    saveFigures(settings.resultPath)

    %HDF5
    nf = input('Would yo create an HDF5 file? y/n: ','s');
    if strcmp(nf,'y')
        buildLuGRE_HDF5([], settings)
    end


    % --- Command Window Log
    fprintf('Complete catalogue elaboration ended. \n\n')
end

end