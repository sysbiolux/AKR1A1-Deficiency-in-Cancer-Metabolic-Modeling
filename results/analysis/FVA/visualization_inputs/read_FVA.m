function [pathways, scores] = read_FVA(filename)

    if ~exist(filename, 'file')
        error('Input file not found: %s', filename);
    end

    % Read Excel file using xlsread for MATLAB R2019b compatibility
    [~, ~, raw] = xlsread(filename);

    if isempty(raw)
        error('The file is empty or could not be read: %s', filename);
    end

    if size(raw, 2) < 2
        error('The file must contain at least two columns: %s', filename);
    end

    % Assume first row contains headers
    header = raw(1, :);

    % Detect pathway column
    pathwayCol = [];

    for j = 1:numel(header)

        headerName = strtrim(char(string(header{j})));

        if strcmpi(headerName, 'Pathways') || strcmpi(headerName, 'subSys')
            pathwayCol = j;
            break;
        end

    end

    % If no recognized header is found, use first column
    if isempty(pathwayCol)
        pathwayCol = 1;
    end

    % Similarity score is assumed to be the second column
    scoreCol = 2;

    % Read data below the header row
    pathwayCells = raw(2:end, pathwayCol);
    scoreCells = raw(2:end, scoreCol);

    pathways = strings(numel(pathwayCells), 1);
    scores = nan(numel(scoreCells), 1);

    for r = 1:numel(pathwayCells)

        % Pathway names
        if isempty(pathwayCells{r})
            pathways(r) = "";
        else
            pathways(r) = string(pathwayCells{r});
        end

        % Numeric scores
        value = scoreCells{r};

        if isnumeric(value)
            scores(r) = value;
        elseif ischar(value) || isstring(value)
            scores(r) = str2double(string(value));
        else
            scores(r) = NaN;
        end

    end

    % Remove empty pathway names and invalid scores
    valid = ...
        ~isnan(scores) & ...
        strlength(strtrim(pathways)) > 0;

    pathways = pathways(valid);
    scores = scores(valid);

    if isempty(scores)
        error('No valid pathway scores were found in: %s', filename);
    end

end