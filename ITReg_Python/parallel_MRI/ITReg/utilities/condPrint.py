
"""
function condPrint(bool,varargin)
if bool
    if isnumeric(varargin{1})
        fprintf(varargin{1},varargin{2:end});
        if varargin{1}>=3
            fprintf(varargin{2:end});
        end
    else
        fprintf(varargin{:})
    end
end
end
"""


def condPrint(b,*varargin):
    '''
    if b:
        if isnumeric(varargin{1})
            fprintf(varargin{1},varargin{2:end});
            if varargin{1}>=3
                fprintf(varargin{2:end});
            end
        else
            fprintf(varargin{:})
        end
    end
    end
    '''
    #I'll just print this stuff. Python should handle it.
    print(varargin)
