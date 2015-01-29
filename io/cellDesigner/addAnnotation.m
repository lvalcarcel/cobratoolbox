function [fname_out,var] = annotation(fname, fname_out,infix,model)


% fname      the original XML file to be modified to include the annotations
% fanme_out  the output file name of the XML file
% infix      the ID (e.g., the metabolite or reaction IDs)
% model      a COBRA model that contains the annotations which can be
%            retrieved by using the infix as the index value.

% e.g., [a,var]=annotation('example_annot.xml','anno_test.xml','r1922',recon2);

% e.g., [a,var]=annotation('example_annot.xml','anno_test.xml',{'r1922','10fthf[c]'},recon2)

% e.g., [a,var]=addAnnotation_b('listOfMets.xml','listOfMets_annoted.xml',{'10fthf5glu[c]'},recon2)


if nargin<4
    
    model=recon2;
end


if nargin<3
    
    prefix='id="';
    infix='r1922';
    suffix='"';
    rxnName=[prefix,infix,suffix];
    
end


prefix='name="';  % for metabolites in recon2


% prefix='id="';

%

suffix='"';

% rxnName={};



for d=1:length(infix)
    
    rxnName(d)=strcat(prefix,infix(d),suffix);
    
end




if nargin<2 || isempty(fname_out)
    
    [fname_out, fpath]=uiputfile('*.xml','CellDesigner SBML Source File');
    if(fname_out==0)
        return;
    end
    f_out=fopen([fpath,fname_out],'w');
else
    f_out=fopen(fname_out,'w');
end

if nargin<1 || isempty(fname)
    [fname, fpath]=uigetfile('*.xml','CellDesigner SBML Source File');
    if(fname==0)
        return;
    end
    f_id=fopen([fpath,fname],'r');
else
    f_id=fopen(fname,'r');
    
end


numOfLine=0;

% rem=fgets(f_id); numOfLine=numOfLine+1;

%%% the template for online map system to recognize

preTxt(1).str='<notes>';
preTxt(2).str='<html xmlns="http://www.w3.org/1999/xhtml">';
preTxt(3).str='<head>';
preTxt(4).str='<title/>';
preTxt(5).str='</head>';
preTxt(6).str='<body>';

preTxt(7).str='</body>';
preTxt(8).str='</html>';
preTxt(9).str='</notes>';

% rem=fgets(f_id);

h = waitbar(0,'Progressing');

MainTxt={};


while ~feof(f_id);
      
    numOfLine=numOfLine+1;
    rem=fgets(f_id);
    %     try
    MainTxt(numOfLine,1)=cellstr(rem);
    
    %     catch
    %         disp(rem);
    %     end
    
end


total_length=length(MainTxt);

disp(total_length);

n=0;  % the line number of the code in the new file that is copied from the original file.
t=1;


for i=1:10;
    ct(i,1)=i*total_length/10;
end

   % met_str='name="'
    
    met_str='<species metaid="'
    
    rxn_str='<reaction metaid="'


for t=1:total_length   % go through each line of the SBML file.
    
    
    if ismember(t, ct)~=0||t==total_length;
        
        disp(t)
        waitbar(t/total_length,h);
        
    end
    
    
    
    n=n+1;
    
    MainTxt_new(n,1)=MainTxt(t);
    

    
    if ((~isempty(strfind(MainTxt(t),met_str)))||(~isempty(strfind(MainTxt(t),met_str))));

        
        for in=1:length(rxnName);% for in=1:length(infix); go though each line of the Rxn List.
            
            % disp(length(rxnName));
            
            line_st=strfind(MainTxt(t),rxnName{in});
            
            
            if ~isempty(line_st{1})  %isempty(line_st)~=0; % the line contains the rxn keywords
                
                
                % msgbox('reaction found');
                disp(line_st)
                
                
                
                %%%%%% preTxt (1:6)
                for p=1:6;
                    
                    MainTxt_new(n+p,1)=cellstr(preTxt(p).str);
                    
                end
                n=n+p;
                total_length=total_length+p;
                
                
                [rxnItems,rxnContent]=contructItems(infix(in),model);
                
                
                for k=1:length(rxnContent);
                    
                    rxnContent(k)          
                    disp(n);
                    disp(k);
                    
                    disp(MainTxt_new(n,1));
                    disp('%%%%%%%%%%%%');
                    disp(rxnItems(k));
                    
                    
                    MainTxt_new(n+k,1)=rxnItems(k);
                    
                    
                    
                end
                
                
                n=n+k;
                total_length=total_length+k;
                
                
                for p_e=1:3
                    
                    MainTxt_new(n+p_e,1)=cellstr(preTxt(p_e+6).str);
                    
                end
                n=n+p_e;
                total_length=total_length+p_e;
                
            end
            
        end
        
    end
    
end


var=MainTxt_new;


for ww=1:length(MainTxt_new);
    
    
    fprintf(f_out,'%s\n',char(MainTxt_new(ww)));
    
    
    
end


close(h);

fclose(f_out);



end


% <body>
% Symbol:
% Abbreviation: ABREVIATION
% Formula: FORMULA
% MechanicalConfidenceScore: MCS
% LowerBound: LB
% UpperBound: UB
% Subsystem: SUBSYSTEM
% GeneProteinReaction: GPR (using miriam registry format, for example: &quot;(urn:miriam:ccds:CCDS8639.1) or (urn:miriam:ccds:CCDS26.1) or (urn:miriam:ccds:CCDS314.2) or (urn:miriam:ccds:CCDS314.1)&quot;)
% Synonyms:
% DESCRIPTION
% </body>



%%%%%%%%%%%%%%%%%%%%%%%%


function [rxnItems,rxnContent]=contructItems(index,model)

% rxnItems    an array of combined Keywords and Contents

% rxnConent   an array of Conetens


infix=index;


%%



num=find(strcmp(model.rxns(:,1),infix)); %% find the reaction number in the model.
para='rxn';


if isempty(num)
    num=find(strcmp(model.mets(:,1),infix));
    para='met';
    
elseif ~isempty(find(strcmp(model.mets(:,1),infix), 1))        %%%%% Error! the reaction and metabolite use the same name.
    
    msg=strcat(infix, ' is used as a reaction name as well as a metabolite name');
    warndlg(msg,'Warning!');
    
end



%% the annotation template for metabolites

if strcmp(para,'met')
    
    
    metKeywords={
        'Symbol: ';
        'Abbreviation: ';
        'ChargedFormula: ';
        'Charge: ';
        'Synonyms: ';
        'Description: '
        };
    
    
    
    %%% assign a initial value of ' ' to the list of the variables.
    
    Symbol=' ';
    Abbreviation=' ';
    ChargedFormula=' ';
    Charge=' ';
    Synonyms=' ';
    Description=' ';
    
    
    
    if ~isfield(model,'Symbol');
        Symbol=' ';
    end
    
    Abbreviation=model.mets(num);
    ChargedFormula=model.metFormulas(num);
    
    
    Charge=model.metCharge(num);
    
    if isnumeric(Charge)
        Charge=num2str(Charge);
    end
    
    
    if ~isfield(model,'Synonyms');
        Synonyms=' ';
    end
    
    
    if ~isfield(model,'Description: ')
        Description=' ';
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    metContent=[
        Symbol(1);
        Abbreviation(1);
        ChargedFormula(1);
        Charge(1);
        Synonyms(1);
        Description(1);
        ];
    
    
    %%%%%
    
    
    for f=1:length(metContent);
        
        rxnItems(f)=strcat(metKeywords(f),metContent(f));
    end
    
    
    rxnContent=metContent;
elseif strcmp(para,'rxn')
    
    
    
    %% check if annotation field for reactions of the COBRA model are avaliable and use a initial value of ' ' as the default annotation.
    
    Abbreviation=model.rxns(num);
    Description=model.rxnNames(num);
    
    if ~isfield(model,'MCS');
        MCS=' ';
    end
    
    
    Ref=' ';
    
    ECNumber=model.rxnECNumbers(num);
    KeggID=model.rxnKeggID(num);
    
    LastModified=' ';
    
    LB=model.lb(num);
    UB=model.ub(num);
    
    if isnumeric(LB)|isnumeric(UB)
        LB=num2str(LB);
        UB=num2str(UB);
        
    end
    
    
    if ~isfield(model,'grRules');
        grRules=' ';
    end
    
    if exist('CS')~=1;
        CS=1;
    end
    
    if isnumeric(CS)
        CS=num2str(CS);
    end
    
    
    GPR=model.grRules(num);
    
    Subsystem=model.subSystems(num);
    
    %%
    
    % MCS-CS: Mechanical Confidence Score CS: Confidence Score - LB: Lower
    % Bound - UB: Upper Bound MCS: Mechanical Confidence Score - GPR: Gene
    % Protein Reaction)
    
    %% the annotation template for metabolites
    
    rxnKeywords={
        'Abbreviation: ';
        'Description: ';
        'MCS: ';
        'Ref: ';
        'ECNumber: ';
        'KeggID: ';
        'LastModified: ';
        'LB: ';
        'UB: ';
        'CS: ';
        'GPR: ';
        'Subsystem: '};
    
    
    
    rxnContent=[
        Abbreviation(1);
        Description(1);
        MCS(1);
        Ref(1);
        ECNumber(1);
        KeggID(1);
        LastModified(1);
        LB;
        UB;
        CS;
        GPR;
        Subsystem];
    
    
    %%
    
    for d=1:length(rxnContent);
        
        rxnItems(d)=strcat(rxnKeywords(d),rxnContent(d));
    end
end
end


