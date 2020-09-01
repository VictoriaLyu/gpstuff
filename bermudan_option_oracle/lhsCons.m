function y = lhsCons(len, rec)

% Generates lhs design in a constrained space

% Description:
%   LHSCONS generates LHS design of length len in a constrained space with 
%   bounds rec.
%   Y = LHSCONS(LEN, REC) returns the constrained LHS design.
%   
%       LEN - length of the generated design
%       REC - constraints of the design

d = size(rec,1);
y = lhsdesign(len,d);

start = repmat(rec(:,1)',len,1);
ed = repmat(rec(:,2)',len,1);
diff = ed - start;
y = y.*diff + start;
end

