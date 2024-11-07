% The COBRAToolbox: testPanModels.m
%
% Purpose:
%     - tests that pan-species models generated by the createPanModels
%     function can produce biomass and produce reasonable amounts of ATP
%
% Author:
%     - Almut Heinken - September 2021
%
% Note:
%     - The solver libraries must be included separately

% initialize the test
fileDir = fileparts(which('testCreatePanModels'));
cd(fileDir);

modelList={
        'Abiotrophia_defectiva_ATCC_49176'
        'Acidaminococcus_fermentans_DSM_20731'
        'Acidaminococcus_intestini_RyC_MR95'
        'Acidaminococcus_sp_D21'
        'Acinetobacter_baumannii_AB0057'
        'Acinetobacter_calcoaceticus_PHEA_2'
        'Acinetobacter_haemolyticus_NIPH_261'
        'Acinetobacter_johnsonii_SH046'
        'Acinetobacter_junii_SH205'
        'Acinetobacter_lwoffii_WJ10621'
        'Acinetobacter_pittii_ANC_4052'
        'Acinetobacter_radioresistens_NIPH_2130'
        };
    for i=1:length(modelList)
        model = getDistributedModel([modelList{i} '.mat']);
        % Save all the models into the modelFolder
        save(fullfile(fileDir, modelList{i}), 'model');
    end
    
numWorkers=4;

% create the pan-models on species level
panPath=[pwd filesep 'panSpeciesModels'];

builtTaxa = {'Abiotrophia defectiva','Acidaminococcus fermentans','Acidaminococcus intestini','Acidaminococcus sp. D21','Acinetobacter_baumannii','Acinetobacter haemolyticus','Acinetobacter johnsonii','Acinetobacter junii','Acinetobacter lwoffii','Acinetobacter pittii','Acinetobacter radioresistens'};
createPanModels(fileDir,panPath,'Species','AGORA_infoFile.xlsx',numWorkers,builtTaxa)

% test that pan-models can grow
[notGrowing,Biomass_fluxes] = plotBiomassTestResults(panPath, 'pan-models','numWorkers',numWorkers);
assert(isempty(notGrowing))

% test that ATP production is not too high
[tooHighATP,ATP_fluxes] = plotATPTestResults(panPath, 'pan-models','numWorkers',numWorkers);
assert(max(cell2mat(ATP_fluxes(2:end,2))) < 200)
assert(max(cell2mat(ATP_fluxes(2:end,3))) < 150)

% create the pan-models on genus level
panPath=[pwd filesep 'panGenusModels'];

builtTaxa = {'Abiotrophia','Acidaminococcus','Acinetobacter'};
createPanModels(fileDir,panPath,'Genus','AGORA_infoFile.xlsx',numWorkers,builtTaxa);

% test that pan-models can grow
[notGrowing,Biomass_fluxes] = plotBiomassTestResults(panPath, 'pan-models','numWorkers',numWorkers);
assert(isempty(notGrowing))

% test that ATP production is not too high
[tooHighATP,ATP_fluxes] = plotATPTestResults(panPath, 'pan-models','numWorkers',numWorkers);
assert(max(cell2mat(ATP_fluxes(2:end,2))) < 250)
assert(max(cell2mat(ATP_fluxes(2:end,3))) < 200)
