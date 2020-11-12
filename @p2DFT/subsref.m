function out = subsref(obj,S)
%function to expose members 
out = builtin('subsref',obj,S);
end

% source:https://stackoverflow.com/questions/64716336/how-to-make-a-member-public-in-an-class-object-created-using-matlab-class-folder/64719282?noredirect=1#comment114455432_64719282

