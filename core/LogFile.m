classdef LogFile < handle
    %UNTITLED7 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        filename
        filepath
        string
        permission
    end
    
    methods
        function obj = LogFile(filepath,varargin)
            IP = inputParser;
            addOptional(IP,'permission','w',@(x)ischar(x))
            parse(IP,varargin{:})
            obj.permission = IP.Results.permission;
            [~,name,~] = fileparts(filepath);
            obj.filename = name;
            obj.filepath = filepath;
        end
        
        function nbytes  = write(obj,str,varargin)
            IP = inputParser;
            addOptional(IP,'permission',obj.permission,@(x)ischar(x))
            parse(IP,varargin{:})
            obj.permission = IP.Results.permission;
            
            fid = fopen(obj.filepath,obj.permission);
            cleanfid = onCleanup(@() fclose(fid));  %called on completion of function, whether normally, by exception or abort
            nbytes = fprintf(fid,str);
            obj.string{end+1} = str;
        end
    end
    
    methods (Static)
        function nbytes = WriteLog(filepath,str,permission)
            fid = fopen(filepath, permission);
            cleanfid = onCleanup(@() fclose(fid));  %called on completion of function, whether normally, by exception or abort
            nbytes = fprintf(fid,str); 
            
            
        end
    end
end

