function [ ae ] = globalTrans( dx,dy )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

ae = [dx dy 0 0 0 0;
      -dy dx 0 0 0 0;
      0 0 1 0 0 0;
      0 0 0 dx dy 0;
      0 0 0 -dy dx 0;
      0 0 0 0 0 1];
  
end

