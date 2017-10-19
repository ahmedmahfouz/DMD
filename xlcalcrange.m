function xlrange = xlcalcrange(refCell,r,c,m,n)
% XLCALCRANGE calculate full target range, in Excel A1 notation.
%   XLRANGE = XLCALCRANGE(REFCELL,R,C,M,N) returns the full target range 
%   XLRANGE in Excel A1 notation to cover an array of size M rows & N columns
%   and starting from an offset of R rows & C columns from REFCELL.
%   Any workbook\worksheet prefix is removed from REFCELL, as are any 
%   absolute '$' markers prior to range calculation. i.e.
%   REFCELL = 'C:\directorypath\[filename.xlsm]worksheet1'!$P$6:$AC$91 
%   becomes:
%   REFCELL = P6:AC91
%
%   XLRANGE = XLCALCRANGE(REFCELL) REFCELL defines an Excel named range, 
%   which is evaluated to obtain the underlying A1 notation; offsets R & C 
%   are set to zero and the array size MxN is set to 1x1
%
%   XLRANGE = XLCALCRANGE(REFCELL,M) REFCELL is a single cell from which
%   XLRANGE becomes an array of size Mx1 with zero offsets.
%
%   XLRANGE = XLCALCRANGE(REFCELL,R) REFCELL is a range of cells and
%   XLRANGE is identically sized, but offset by R rows.
%
%   XLRANGE = XLCALCRANGE(REFCELL,M,N) REFCELL is a single cell from which
%   XLRANGE becomes an array of size MxN with zero offsets.
%
%   XLRANGE = XLCALCRANGE(REFCELL,R,C) REFCELL is a range of cells and
%   XLRANGE is identically sized, but offset by R rows and C columns.
%
%   XLRANGE = XLCALCRANGE(REFCELL,R,C,M) REFCELL can be either a single cell
%   or a cell range.  The output XLRANGE has the same column width as
%   REFCELL, but has row height M and is offset from REFCELL by R rows and 
%   C columns.
%
%   XLRANGE = XLCALCRANGE(REFCELL,R,C,M,N)  REFCELL can be either a single 
%   cell or a cell range.  The output XLRANGE has row height M and column
%   width N and is offset from REFCELL by R rows and C columns.
%
%
%   INPUT PARAMETERS:
%      refCell:  string defining the top left-hand start of the range to
%                be calculated
%                string defining an Excel named range
%      r:        double defining a whole number of offsetting rows to add/
%                subtract to/from the refCell row number
%      c:        double defining a whole number of offsetting columns to
%                add/subtract to/from the refCell column
%      m:        double defining the whole number of rows in the range to
%                be calculated
%      n:        double defining the whole number of columns in the range
%                to be calculated
%
%   RETURN PARAMETERS:
%      xlrange:  string containing the full target range in Excel A1
%                notation
%
%   by Richard de Garis
%   Version 1.0   08/02/2010 - first version
%   Version 1.0.1 08/10/2010 - fixed bug that requires removing workbook/
%                              worksheet reference from rangeLabel.  This
%                              is important when named ranges are
%                              duplicated when worksheets are copied.
%   Version 2.0   08/20/2010 - allowed more flexibility in working with an
%                              input cell or cell range, which drives the
%                              interpretation of subsequent input values.
%                              Also makes use of a custom function that
%                              checks the validity of an Excel range.
%   Version 2.0.1 08/30/2010 - changed function name from calcxlrange to
%                              xlcalcrange to be consistent with all other
%                              functions developed for use with Excel.
%   Version 2.0.2 12/13/2010 - corrected a bug in the way offsets are
%                              applied to the LRH row-col pair, if it
%                              exists - the offsets were applied twice.
%                              Now the LRH row-col pair simply reflects the
%                              ULH row-col pair.
%==========================================================================

if ~ischar(refCell)
    error('MATLAB:RdeGFunction:xlcalcrange:nonStringInput',...
        'Excel range must be a character string');
end

% initialize the output
xlrange = '';

RC = xlvalidaterange(refCell);

if isequal(RC,0)
    % input is potentially an Excel named range
    xlrange = xlgetnamedrange(refCell);
    if isequal(xlrange,'')
        error('MATLAB:RdeGFunction:xlcalcrange:invalidInput',...
            'Excel range must either be explicit or a valid named range');
    else
        RC = xlvalidaterange(xlrange);
    end
end

nRC = size(RC,1);
switch nRC
    case 1
        cellRange = false;
    case 2
        if isequal([RC{1,1} RC{1,2}],[RC{2,1} RC{2,2}])
            cellRange = false;
        else
            cellRange = true;
        end
end

if nargin < 2
    r = 0;
    c = 0;
    m = 0;
    n = 0;
    if nRC >= 1
        xlrange = [RC{1,1} num2str(RC{1,2})];
    end
    if cellRange
        xlrange = [xlrange ':' RC{2,1} num2str(RC{2,2})];
    end
else if nargin < 3
        if cellRange
            % cell range with row-offset, r, specified
            c = 0;
            m = 0;
            n = 0;
        else
            % single cell with output row-size, m, specified
            m = r;
            n = 0;
            r = 0;
            c = 0;
        end
    else if nargin < 4 
            if cellRange
                % cell range with offsets, r, c, specified
                m = 0;
                n = 0;
            else
               % single cell with output size, m, n, specified
                m = r;
                n = c;
                r = 0;
                c = 0;
            end
        else if nargin < 5
                % Both offets, r, c & row-size, m, specified
                n = 0;
            end
        end
    end
end

if m < 0 || n < 0
    error('MATLAB:RdeGFunction:xlcalcrange:invalidArraySize',...
        'Array dimensions MxN must be greater than or equal to zero');
end

if m ~= floor(m) || n ~= floor(n) || r ~= floor(r) || c ~= floor(c)
    error('MATLAB:RdeGFunction:xlcalcrange:nonInteger',...
        'Non-integer values are not allowed');
end

try
    % Apply the offsets
    RC{1,1} = dec2base27(base27dec(RC{1,1}) + c);
    RC{1,2} = RC{1,2} + r;
    if nRC > 1
        RC{2,1} = dec2base27(base27dec(RC{2,1}) + c);
        RC{2,2} = RC{2,2} + r;
    else
        RC{2,1} = RC{1,1};
        RC{2,2} = RC{1,2};
    end
    
    % Apply the array sizing
    if n > 0
        RC{2,1} = dec2base27(base27dec(RC{1,1}) + n - 1);
    end
    
    if m > 0
        RC{2,2} = RC{1,2} + m - 1;
    end
    
    % Check validation of resulting range
    xlrange = [RC{1,1} num2str(RC{1,2}) ':' RC{2,1} num2str(RC{2,2})];
    RC      = xlvalidaterange(xlrange);
    if isequal(RC,0)
        error;
    end
catch exception
    error('MATLAB:RdeGFunction:xlcalcrange:invalidNewRange',...
        'Offsets and/or array size are invalid for the input Excel range.');
end
end

%--------------------------------------------------------------------------


function xlrange = xlgetnamedrange(name)

xlrange = '';
try
    Excel      = evalin('base', 'Excel');
    xlNameList = Excel.ActiveWorkbook.Names;
    nNames     = xlNameList.count;

    for i = 1:nNames
        xlNamedRange = xlNameList.Item(i);
        xlRangeLabel = xlNamedRange.Name;
        % remove workbook/worksheet prefix
        xlRangeLabel(1:strfind(xlRangeLabel,'!')) = []; 
        if(strcmpi(xlRangeLabel,name)) 
            xlrange = get(xlNamedRange,'RefersToLocal');
            break
        end
    end
    error;
catch exception
    if ~exist('Excel', 'var')
        error('MATLAB:RdeGFunction:xlcalcrange:inactiveExcelServer',...
            'Excel COM Server must be running to evaluate named range');
    end
    if ~exist('xlNameList','var')
        error ('MATLAB:RdeGFunction:xlcalcrange:noOpenWorkbook',...
            'There must be an open Excel workbook to evaluate named range');
    end
    if isequal(xlrange, '')
        error('MATLAB:RdeGFunction:xlcalcrange:invalidNamedRange',...
            'Specified named range must be valid for the active Excel workbook');
    end
end
end




function rc = xlvalidaterange(strIn)
% XLVALIDATERANGE validate legality of Excel range that is in A1 notation.
%   RC = XLVALIDATERANGE(STRIN) given STRIN is a valid Excel range, for 
%   Excel versions 2007 and up, returns a cell array of size n x 2 
%   containing the disaggregated row and column components, in A1
%   notation.  n equals 1 for a single cell and n equals 2 for a cell range.
%   Returns zero if STRIN is an invalid Excel range.  Excel named ranges
%   are treated as invalid.
%   Any workbook\worksheet prefix is removed from STRIN, as are any 
%   absolute '$' markers prior to validation.
%
%    INPUT PARAMETERS:
%       strIn:  string describing an Excel cell e.g. 'A1', or cell range
%               e.g. 'A1:B2'.
%
%    RETURN PARAMETERS:
%       rc:     double containing zero in the case of an invalid Excel
%               range.
%       rc:     cell array of size (n,2) containing separated row and 
%               column components of valid Excel cell (n=1) or cell range
%               (n=2). e.g. {'A', 1; 'B', 2}.  The column component is a
%               string.  The row component is a double.
%
%   by Richard de Garis
%   Version 1.0   08/20/2010 - first version
%==========================================================================

%% Set up constants
MIN_RC     = 1;
MAX_ROW    = 1048576;
MAX_COL    = base27dec('XFD');
RC_PATTERN = '\$?[A-Z]{1,3}\$?\d{1,7}';

% End Set up constants

%% Error chech the input

if ~ischar(strIn)
    error('MATLAB:RdeGFunction:xlvalidaterange:nonStringInput',...
        'Excel range must be a character string');
else
    strIn(1:strfind(strIn,'!')) = []; % remove workbook/worksheet prefix
    strIn                       = upper(strIn);
end

% End Error chech the input

%% Parse input
RCmatch = regexp(strIn,RC_PATTERN,'match');
nMatch  = length(RCmatch);
nRC     = 0;
rc      = 0;

switch nMatch
    case 1
        if isequal(strIn, RCmatch{1})
            % input string is potentially a valid single cell
            nRC = 1;
        end
    case 2
        if isequal(strIn, [RCmatch{1} ':' RCmatch{2}])
            % input string is potentially a valid cell range
            nRC = 2;
        end
end

if nRC > 0
    % Check input for legal Excel Range
    RCmatch = strrep(RCmatch,'$','');
    rc      = cell(nRC,2);
    for i = 1:nRC
        colA = RCmatch{i}(isletter(RCmatch{i}));
        colN = base27dec(colA);
        row  = str2double(RCmatch{i}(~isletter(RCmatch{i})));

        if colN <= MAX_COL && colN >= MIN_RC &&...
            row <= MAX_ROW && row >= MIN_RC
            rc{i,1} = colA;
            rc{i,2} = row;
        else
            rc = 0;
            break;
        end
    end
end


% End Parse input
end




function s = dec2base27(d)

%   DEC2BASE27(D) returns the representation of D as a string in base 27,
%   expressed as 'A'..'Z', 'AA','AB'...'AZ', and so on. Note, there is no zero
%   digit, so strictly we have hybrid base26, base27 number system.  D must be a
%   negative integer bigger than 0 and smaller than 2^52.
%
%   Examples
%       dec2base(1) returns 'A'
%       dec2base(26) returns 'Z'
%       dec2base(27) returns 'AA'
%
%   Version 2.0
%==========================================================================
%   Based on the function DEC2BASE27 contained within the Matlab function
%   XLSWRITE.
%   Copyright 1984-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.17 $  $Date: 2009/11/16 22:26:56 $
%==========================================================================

d = d(:);
if d ~= floor(d) || any(d(:) < 0) || any(d(:) > 1/eps)
    error('MATLAB:xlswrite:Dec2BaseInput',...
        'D must be an integer, 0 <= D <= 2^52.');
end

n_letters   = 1;
begin       = 0;
current_sum = 26;
letters     = 'A':'Z';

while d > current_sum
    n_letters   = n_letters + 1;
    begin       = current_sum;
    current_sum = begin + 26.^n_letters;
end

linear_index = d - begin;
outputs      = cell(1,n_letters);

% Here, ind2sub finds the equivalent subscripts of the linear index in an 
% array of predefined dimensions in multiplexs of 26 i.e. 26, 26x26,
% 26x26x26, etc.
% For example the linear index of an array of 26x26 goes from 1 in the top
% left-hand corner, increases sequentially down the rows, top-to-bottom,
% and then continues along the columns, left-to-right.  Thus the subscripts
% of:
% 1 are [1,1], 
% 26 are [26,1],
% 27 are [1,2],
% 676 are [26,26]
[outputs{1:n_letters}] = ind2sub(repmat(26,1,n_letters), linear_index);
s                      = fliplr(letters([outputs{:}]));
end





function d = base27dec(s)
%   BASE27DEC(S) returns the decimal of string S which represents a number in
%   base 27, expressed as 'A'..'Z', 'AA','AB'...'AZ', and so on. Note, there is
%   no zero so strictly we have hybrid base26, base27 number system.
%
%   Examples
%       base27dec('A') returns 1
%       base27dec('Z') returns 26
%       base27dec('IV') returns 256
%
%   Version 2.0
%==========================================================================
%   Based on the function BASE27DEC contained within the Matlab function
%   XLSWRITE.
%   Copyright 1984-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.17 $  $Date: 2009/11/16 22:26:56 $
%==========================================================================

s = upper(s);
if length(s) == 1
   d = s(1) -'A' + 1;
else
    cumulative = 0;
    for i = 1:numel(s)-1
        cumulative = cumulative + 26.^i;
    end
    indexes_flipped = 1 + s - 'A';
    indexes = fliplr(indexes_flipped);
    indexes_in_cells = mat2cell(indexes, 1, ones(1,numel(indexes)));
    d = cumulative + sub2ind(repmat(26, 1,numel(s)), indexes_in_cells{:});
end
end
