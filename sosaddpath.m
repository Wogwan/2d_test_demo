% Add paths to SOS analysis toolboxes
%

% Add multipoly

cm = computer;

if cm(1) == 'M' ||  cm(1)=='G'
    set(0,'DefaultFigureWindowStyle','docked')
    % Add chebfun-master
    addpath([pwd '/toolbox/chebfun-master']);
    
    % Add gpml-matlab-master
    addpath([pwd '/toolbox/gpml']);
    run([pwd '/toolbox/gpml/startup.m'])
    
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
    
    % Add utils
    addpath([pwd '/utils'])
    
    % Add updates
    addpath([pwd '/updates'])
    
    % Add demo_2D
    addpath([pwd '/demo_2D'])
    addpath([pwd '/demo_2D/paper_example'])
    
    % Add demo_3D
    addpath([pwd '/demo_3D'])
    addpath([pwd '/demo_3D/paper_example'])
    
elseif cm(1) == 'P'
    set(0,'DefaultFigureWindowStyle','docked')
%     set(0,'DefaultFigureWindowStyle','normal')
    % Add chebfun-master
    addpath([pwd '/toolbox/chebfun-master']);
    
    % Add gpml-matlab-master
    addpath([pwd '/toolbox/gpml']);
    run([pwd '/toolbox/gpml/startup.m'])
    
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
    
    % Add utils
    addpath([pwd '/utils'])
    
    % Add updates
    addpath([pwd '/updates'])
    
    % Add demo_2D
    addpath([pwd '/demo_2D'])
    addpath([pwd '/demo_2D/paper_example'])
    
    % Add demo_3D
    addpath([pwd '/demo_3D'])
    addpath([pwd '/demo_3D/paper_example'])
end


