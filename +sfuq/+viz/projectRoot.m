function p = projectRoot()
%PROJECTROOT  Absolute path to the repository root.
%
%   Derived from this file's own location, so every script works regardless of
%   the current working directory. Hard-coding relative paths is the most
%   common reason a cloned research repository fails on someone else's
%   machine.

    here = fileparts(mfilename('fullpath'));      % .../+sfuq/+viz
    p = fileparts(fileparts(here));               % repository root
end
