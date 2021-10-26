% Add paths to SOS analysis toolboxes
%

% Add multipoly

cm = computer;

if cm(1) == 'M' ||  cm(1)=='G'
%     % Add chebfun-master
%     addpath([pwd '/toolbox/chebfun-master']);
    
    % Add gpml-matlab-master
    addpath([pwd '/toolbox/gpml']); 
    % Please run up the '/gpml-matlab-master/startup.m' if you first
    % running the code
    
    % Add multipoly
    addpath([pwd '/toolbox/multipoly']);
    
    % Add nlanal
    addpath([pwd '/toolbox/nlanal']);

    % Add my version of SOSTools
    addpath([pwd '/toolbox/sosopt']);
    addpath([pwd '/toolbox/sosopt/Demos']);

    % Add polysys
    addpath([pwd '/toolbox/polysystems_1_0_3'])
    addpath([pwd '/toolbox/polysystems_1_0_3/demo'])

%     % Add SeDuMi_1_3
%     addpath([pwd '/toolbox/SeDuMi_1_3'])
    
    % Add utils
    addpath([pwd '/utils'])
    
    % Add updates
    addpath([pwd '/update'])
        
    % Add test_function
    addpath([pwd '/test_function'])

elseif cm(1) == 'P'

%     % Add chebfun-master
%     addpath([pwd '/toolbox/chebfun-master']);
    
    % Add gpml-matlab-master
    addpath([pwd '/toolbox/gpml']);
    % Please run up the '/gpml-matlab-master/startup.m' if you first
    % running the code
    
    % Add multipoly
    addpath([pwd '/toolbox/multipoly']);
    
    % Add nlanal
    addpath([pwd '/toolbox/nlanal']);

    % Add my version of SOSTools
    addpath([pwd '/toolbox/sosopt']);
    addpath([pwd '/toolbox/sosopt/Demos']);

    % Add polysys
    addpath([pwd '/toolbox/polysystems_1_0_3'])
    addpath([pwd '/toolbox/polysystems_1_0_3/demo'])

    % Add SeDuMi_1_3
    addpath([pwd '/toolbox/SeDuMi_1_3'])
    
    % Add utils
    addpath([pwd '/utils'])
            
%     Add updates
    addpath([pwd '/update'])
            
    % Add test_function
    addpath([pwd '/test_function'])

end
    

