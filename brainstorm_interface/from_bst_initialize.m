function from_bst_initialize(PARAMS)
    
    % Start Brainstrom if not yet started
    if ~brainstorm('status')
        brainstorm nogui
    end

    % Set the right protocol as current if it is not yet
    protocol_name = PARAMS.protocol_name;
    iProtocol = bst_get('Protocol', protocol_name);
    
    if isempty(iProtocol)
        disp(sprintf('No protocol named ''%s'' could be found', protocol_name));
    elseif bst_get('iProtocol') ~= iProtocol
        gui_brainstorm('SetCurrentProtocol', iProtocol);
        disp(sprintf('Protocol ''%s'' has been activated', protocol_name));
    end

    
end